#!/usr/bin/env python3
"""
Generate country boundary data for the B3D globe example.

Two-stage workflow:
  1. Download: Fetch Natural Earth data and create binary asset
     python3 scripts/gen-globe-data.py download -o assets/globe-countries.bin

  2. Header: Convert binary asset to C header for embedding
     python3 scripts/gen-globe-data.py header -i assets/globe-countries.bin -o examples/globe-data.h

Binary format (delta-encoded with varint compression):
  - uint32: polygon count
  - Per polygon:
    - uint16: point count
    - int16, int16: first point (quantized lon, lat)
    - Subsequent: varint-encoded deltas (ZigZag + LEB128)

Dependencies: pip3 install pyshp requests
"""

from __future__ import annotations

import argparse
import math
import os
import struct
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Optional

try:
    import shapefile  # type: ignore
except ImportError:
    shapefile = None

try:
    import requests  # type: ignore
except ImportError:
    requests = None

# Type aliases for clarity
Point = tuple[float, float]
Polygon = list[Point]

# Natural Earth 1:10m countries shapefile URL
NE_URL = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"

# Default paths
DEFAULT_BINARY_OUTPUT = "assets/globe-countries.bin"
DEFAULT_HEADER_OUTPUT = "examples/globe-data.h"
DEFAULT_TOLERANCE = 0.005
DEFAULT_MAX_POINTS = 4000  # Limit for C source size

# Quantization scale matching globe.c
QUANTIZATION_SCALE = 32767.0


def download_shapefile(url: str, cache_dir: Optional[str] = None) -> str:
    """Download and extract Natural Earth shapefile.

    @url: URL to download from
    @cache_dir: directory for caching (default: system temp)
    Returns path to extracted .shp file.
    Raises SystemExit on network or extraction errors.
    """
    if not requests:
        print("Error: requests module required for downloading.")
        print("Install with: pip3 install requests")
        sys.exit(1)

    if cache_dir is None:
        cache_dir = tempfile.gettempdir()

    zip_path = os.path.join(cache_dir, "ne_10m_countries.zip")
    extract_dir = os.path.join(cache_dir, "ne_10m_countries")

    shp_path = os.path.join(extract_dir, "ne_10m_admin_0_countries.shp")
    if os.path.exists(shp_path):
        print(f"Using cached shapefile: {shp_path}")
        return shp_path

    print(f"Downloading Natural Earth data from {url}...")
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error downloading data: {e}")
        sys.exit(1)

    try:
        with open(zip_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    except IOError as e:
        print(f"Error writing zip file: {e}")
        sys.exit(1)

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Verify zip integrity before extraction
            bad_file = zf.testzip()
            if bad_file:
                print(f"Error: Corrupted file in zip: {bad_file}")
                os.remove(zip_path)  # Remove corrupted zip for re-download
                sys.exit(1)
            zf.extractall(extract_dir)
    except zipfile.BadZipFile as e:
        print(f"Error extracting zip file: {e}")
        if os.path.exists(zip_path):
            os.remove(zip_path)  # Remove corrupted zip for re-download
        sys.exit(1)

    if not os.path.exists(shp_path):
        print(f"Error: Expected shapefile not found: {shp_path}")
        sys.exit(1)

    return shp_path


def douglas_peucker(points: list[Point], tolerance: float) -> list[Point]:
    """Simplify a polyline using the Douglas-Peucker algorithm.

    @points: list of (x, y) tuples
    @tolerance: maximum perpendicular distance threshold
    Returns simplified list of points.

    Handles closed polygons (first == last) by simplifying without the
    duplicate endpoint, then restoring it.
    """
    if len(points) <= 2:
        return points

    # Handle closed polygons: remove duplicate endpoint, simplify, restore
    start, end = points[0], points[-1]
    is_closed = start[0] == end[0] and start[1] == end[1]
    if is_closed and len(points) > 3:
        simplified = douglas_peucker(points[:-1], tolerance)
        return simplified + [simplified[0]]

    x1, y1 = start
    x2, y2 = end
    dx, dy = x2 - x1, y2 - y1

    max_dist = 0.0
    max_idx = 0

    denom = math.sqrt(dx * dx + dy * dy) if (dx != 0 or dy != 0) else 1.0

    for i in range(1, len(points) - 1):
        px, py = points[i]
        dist = abs(dy * px - dx * py + x2 * y1 - y2 * x1) / denom
        if dist > max_dist:
            max_dist = dist
            max_idx = i

    if max_dist > tolerance:
        left = douglas_peucker(points[: max_idx + 1], tolerance)
        right = douglas_peucker(points[max_idx:], tolerance)
        return left[:-1] + right
    return [start, end]


def extract_polygons(shp_path: str, tolerance: float = 0) -> list[Polygon]:
    """Extract polygons from shapefile, optionally simplifying them.

    @shp_path: path to .shp file
    @tolerance: Douglas-Peucker simplification tolerance (0 = no simplification)
    Returns list of polygons in radians.
    Raises SystemExit on shapefile read errors.
    """
    if not shapefile:
        print("Error: shapefile module required.")
        print("Install with: pip3 install pyshp")
        sys.exit(1)

    try:
        sf = shapefile.Reader(shp_path)
    except Exception as e:
        print(f"Error reading shapefile: {e}")
        sys.exit(1)

    polygons: list[Polygon] = []

    try:
        for shape in sf.shapes():
            if shape.shapeType not in (
                shapefile.POLYGON,
                shapefile.POLYGONZ,
                shapefile.POLYGONM,
            ):
                continue

            parts = list(shape.parts) + [len(shape.points)]
            for i in range(len(parts) - 1):
                points = shape.points[parts[i] : parts[i + 1]]
                if tolerance > 0 and len(points) > 2:
                    points = douglas_peucker(points, tolerance)

                poly_radians: Polygon = [
                    (math.radians(lon), math.radians(lat)) for lon, lat in points
                ]

                if len(poly_radians) >= 3:
                    polygons.append(poly_radians)
    finally:
        sf.close()

    return polygons


def zigzag_encode(n: int) -> int:
    """Encode signed integer using ZigZag encoding."""
    return (n << 1) ^ (n >> 31)


def write_varint(n: int) -> bytearray:
    """Encode unsigned integer as LEB128 variable-length bytes."""
    out = bytearray()
    while True:
        byte = n & 0x7F
        n >>= 7
        if n == 0:
            out.append(byte)
            break
        out.append(byte | 0x80)
    return out


def poly_area(poly: Polygon) -> float:
    """Calculate polygon area using shoelace formula (for sorting by size)."""
    area = 0.0
    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        area += poly[i][0] * poly[j][1]
        area -= poly[j][0] * poly[i][1]
    return abs(area) / 2.0


def encode_polygons(polygons: list[Polygon]) -> bytearray:
    """Encode polygons to compact binary format with delta encoding.

    Uses rounding for quantization to minimize directional bias.
    Polygons with >65535 points are skipped (uint16 limit).
    """
    out_bytes = bytearray()

    # Filter out oversized polygons
    valid_polys = [p for p in polygons if len(p) <= 65535]
    if len(valid_polys) < len(polygons):
        skipped = len(polygons) - len(valid_polys)
        print(f"Warning: Skipped {skipped} polygons exceeding 65535 points")

    out_bytes.extend(struct.pack("<I", len(valid_polys)))

    for poly in valid_polys:
        out_bytes.extend(struct.pack("<H", len(poly)))
        last_lon_q, last_lat_q = 0, 0

        for i, (lon, lat) in enumerate(poly):
            # Use rounding instead of truncation for better accuracy
            lon_q = int(round(lon / math.pi * QUANTIZATION_SCALE))
            lat_q = int(round(lat / (math.pi / 2.0) * QUANTIZATION_SCALE))
            lon_q = max(-32767, min(32767, lon_q))
            lat_q = max(-32767, min(32767, lat_q))

            if i == 0:
                out_bytes.extend(struct.pack("<hh", lon_q, lat_q))
            else:
                d_lon = lon_q - last_lon_q
                d_lat = lat_q - last_lat_q
                out_bytes.extend(write_varint(zigzag_encode(d_lon)))
                out_bytes.extend(write_varint(zigzag_encode(d_lat)))
            last_lon_q = lon_q
            last_lat_q = lat_q
    return out_bytes


def decode_binary(data: bytes) -> tuple[int, int]:
    """Decode binary header to get polygon and point counts.

    Returns (polygon_count, total_points). Returns (0, 0) on truncated/invalid data.
    """
    if len(data) < 4:
        return 0, 0

    polygon_count = struct.unpack("<I", data[:4])[0]
    offset = 4
    total_points = 0

    # Sanity check: prevent excessive iteration on malformed data
    if polygon_count > 100000:
        return 0, 0

    for _ in range(polygon_count):
        if offset + 2 > len(data):
            return 0, 0  # Truncated data
        point_count = struct.unpack("<H", data[offset : offset + 2])[0]
        total_points += point_count
        offset += 2

        # Skip first point (4 bytes)
        if point_count > 0:
            if offset + 4 > len(data):
                return 0, 0  # Truncated data
            offset += 4

        # Skip varints for remaining points
        for _ in range(1, point_count):
            # Skip lon delta varint
            while offset < len(data) and (data[offset] & 0x80):
                offset += 1
            if offset >= len(data):
                return 0, 0  # Truncated varint
            offset += 1
            # Skip lat delta varint
            while offset < len(data) and (data[offset] & 0x80):
                offset += 1
            if offset >= len(data):
                return 0, 0  # Truncated varint
            offset += 1

    return polygon_count, total_points


def cmd_download(args: argparse.Namespace) -> None:
    """Download Natural Earth data and create binary asset."""
    if not shapefile or not requests:
        print("Required modules (pyshp, requests) not found.")
        print("Install with: pip3 install pyshp requests")
        sys.exit(1)

    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    shp_path = download_shapefile(NE_URL, args.cache_dir)
    polygons = extract_polygons(shp_path, args.simplify)

    print(f"Extracted {len(polygons)} polygons")

    data = encode_polygons(polygons)

    with open(args.output, "wb") as f:
        f.write(data)

    polygon_count, total_points = decode_binary(data)
    print(
        f"Written {len(data)} bytes ({polygon_count} polygons, "
        f"{total_points} points) to {args.output}"
    )


def cmd_header(args: argparse.Namespace) -> None:
    """Convert binary asset to C header."""
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        print("Run 'make update-globe-data' first to download the data.")
        sys.exit(1)

    with open(args.input, "rb") as f:
        data = bytearray(f.read())

    polygon_count, total_points = decode_binary(bytes(data))

    # Apply point limit if specified
    if args.max_points > 0 and total_points > args.max_points:
        print(f"Note: Binary has {total_points} points, limit is {args.max_points}")
        print(
            "Consider regenerating with 'make update-globe-data' using higher simplification"
        )

    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write("/* Generated by scripts/gen-globe-data.py */\n")
        f.write("#include <stdint.h>\n\n")
        f.write(
            f"/* {polygon_count} polygons, {total_points} points, "
            f"size: {len(data)} bytes */\n"
        )
        f.write("static const uint8_t globe_data[] = {\n")
        line = "    "
        for b in data:
            line += f"0x{b:02x},"
            if len(line) > 75:
                f.write(line + "\n")
                line = "    "
        if line.strip():
            f.write(line + "\n")
        f.write("};\n")

    print(
        f"Written C header ({polygon_count} polygons, {len(data)} bytes) "
        f"to {args.output}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate country boundary data for B3D globe example.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download Natural Earth data and create binary asset
  %(prog)s download -o assets/globe-countries.bin

  # Convert binary to C header
  %(prog)s header -i assets/globe-countries.bin -o examples/globe-data.h

  # One-step generation (legacy mode)
  %(prog)s download --format c -o examples/globe-data.h
""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Download subcommand
    dl_parser = subparsers.add_parser(
        "download", help="Download Natural Earth data and create binary asset"
    )
    dl_parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_BINARY_OUTPUT,
        help=f"Output binary file (default: {DEFAULT_BINARY_OUTPUT})",
    )
    dl_parser.add_argument(
        "-s",
        "--simplify",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=f"Simplification tolerance (default: {DEFAULT_TOLERANCE})",
    )
    dl_parser.add_argument(
        "--cache-dir", default=None, help="Directory to cache downloaded shapefiles"
    )
    dl_parser.add_argument(
        "--format",
        choices=["binary", "c"],
        default="binary",
        help="Output format (default: binary)",
    )
    dl_parser.add_argument(
        "--max-points",
        type=int,
        default=DEFAULT_MAX_POINTS,
        help=f"Max points for C format (default: {DEFAULT_MAX_POINTS})",
    )

    # Header subcommand
    hdr_parser = subparsers.add_parser(
        "header", help="Convert binary asset to C header"
    )
    hdr_parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_BINARY_OUTPUT,
        help=f"Input binary file (default: {DEFAULT_BINARY_OUTPUT})",
    )
    hdr_parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_HEADER_OUTPUT,
        help=f"Output C header (default: {DEFAULT_HEADER_OUTPUT})",
    )
    hdr_parser.add_argument(
        "--max-points",
        type=int,
        default=0,
        help="Warn if binary exceeds this point count (0 = no check)",
    )

    args = parser.parse_args()

    if args.command == "download":
        if args.format == "c":
            # Legacy mode: download and generate C header directly
            if not shapefile or not requests:
                print("Required modules (pyshp, requests) not found.")
                print("Install with: pip3 install pyshp requests")
                sys.exit(1)

            Path(args.output).parent.mkdir(parents=True, exist_ok=True)
            shp_path = download_shapefile(NE_URL, args.cache_dir)
            polygons = extract_polygons(shp_path, args.simplify)

            # Filter for embedding
            polygons_with_area = [(poly_area(p), p) for p in polygons]
            polygons_with_area.sort(key=lambda x: x[0], reverse=True)

            selected_polygons: list[Polygon] = []
            total_points = 0

            for _, p in polygons_with_area:
                if total_points + len(p) > args.max_points:
                    continue
                selected_polygons.append(p)
                total_points += len(p)

            print(
                f"Selected {len(selected_polygons)} polygons "
                f"({total_points} points) for embedding."
            )

            data = encode_polygons(selected_polygons)

            with open(args.output, "w", encoding="utf-8") as f:
                f.write("/* Generated by scripts/gen-globe-data.py */\n")
                f.write("#include <stdint.h>\n\n")
                f.write(
                    f"/* {len(selected_polygons)} polygons, "
                    f"size: {len(data)} bytes */\n"
                )
                f.write("static const uint8_t globe_data[] = {\n")
                line = "    "
                for b in data:
                    line += f"0x{b:02x},"
                    if len(line) > 75:
                        f.write(line + "\n")
                        line = "    "
                if line.strip():
                    f.write(line + "\n")
                f.write("};\n")

            print(f"Written to {args.output}")
        else:
            cmd_download(args)
    elif args.command == "header":
        cmd_header(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
