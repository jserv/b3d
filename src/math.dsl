# B3D Math DSL - I❤LA-Style Linear Algebra
#
# Mathematical notation that compiles to C and LaTeX.
# Supports both Unicode symbols and LaTeX-like ASCII escapes.
#
# Syntax:
#   name(params) = expression
#   where
#     param ∈ Type        # or: param \in Type
#
# Notation equivalents (Unicode / ASCII):
#   Types:    ℝ / \R      ℤ / \Z      ⁴ / ^4      ⁴ˣ⁴ / ^4x4
#   Ops:      ∑ / \sum    ∏ / \prod   ∈ / \in
#             · / \cdot   × / \times  √ / \sqrt
#             ‖v‖ / ||v||
#   Greek:    θ / \theta  π / \pi     ε / \eps    δ / \delta
#
# Summation:     ∑(i∈xyz) expr    or  \sum(i \in xyz) expr
# Comprehension: [expr | i∈xyzw] or  [expr | i \in xyzw]

# ═══════════════════════════════════════════════════════════════
# Vector Operations
# ═══════════════════════════════════════════════════════════════

# Dot product: a · b = Σᵢ aᵢbᵢ
vec_dot(a, b) = ∑(i∈xyz) a[i] * b[i]
where
    a ∈ ℝ⁴
    b ∈ ℝ⁴

# Squared length: ‖v‖²
vec_length_sq(v) = vec_dot(v, v)
where
    v ∈ ℝ⁴

# Length: ‖v‖
vec_length(v) = √(vec_dot(v, v))
where
    v ∈ ℝ⁴

# Vector addition: a + b
vec_add(a, b) = [a[i] + b[i] | i∈xyzw]
where
    a ∈ ℝ⁴
    b ∈ ℝ⁴

# Vector subtraction: a - b
vec_sub(a, b) = [a[i] - b[i] | i∈xyzw]
where
    a ∈ ℝ⁴
    b ∈ ℝ⁴

# Scalar multiply: v * s
vec_mul(v, s) = [v[i] * s | i∈xyzw]
where
    v ∈ ℝ⁴
    s ∈ ℝ

# Perspective divide (w=1)
vec_div(v, s) = if |s| < ε
    then [v.x, v.y, v.z, 1]
    else [v.x / s, v.y / s, v.z / s, 1]
where
    v ∈ ℝ⁴
    s ∈ ℝ

# Negate xyz, preserve w
vec_neg(v) = [-v.x, -v.y, -v.z, v.w]
where
    v ∈ ℝ⁴

# Cross product: a × b
vec_cross(a, b) = [
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
    1
]
where
    a ∈ ℝ⁴
    b ∈ ℝ⁴

# Normalize: v / ‖v‖
vec_norm(v) = let len = vec_length(v) in
    if len < ε
    then [0, 0, 0, 1]
    else vec_div(v, len)
where
    v ∈ ℝ⁴

# Linear interpolation: a + t(b - a)
vec_lerp(a, b, t) = [a[i] + t * (b[i] - a[i]) | i∈xyzw]
where
    a ∈ ℝ⁴
    b ∈ ℝ⁴
    t ∈ ℝ

# ═══════════════════════════════════════════════════════════════
# Matrix Construction
# ═══════════════════════════════════════════════════════════════

mat_ident() = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
]

# Extract row as direction (w=0)
mat_row3(M, r) = [M[r][0], M[r][1], M[r][2], 0]
where
    M ∈ ℝ⁴ˣ⁴
    r ∈ ℤ

# Construct matrix from 4 row vectors
mat_from_rows(r0, r1, r2, r3) = [
    [r0.x, r0.y, r0.z, r0.w],
    [r1.x, r1.y, r1.z, r1.w],
    [r2.x, r2.y, r2.z, r2.w],
    [r3.x, r3.y, r3.z, r3.w]
]
where
    r0 ∈ ℝ⁴
    r1 ∈ ℝ⁴
    r2 ∈ ℝ⁴
    r3 ∈ ℝ⁴

# ═══════════════════════════════════════════════════════════════
# Transform Matrices
# ═══════════════════════════════════════════════════════════════

# Rotation about X-axis
mat_rot_x(θ) = [
    [1, 0,       0,      0],
    [0, cos(θ),  sin(θ), 0],
    [0, -sin(θ), cos(θ), 0],
    [0, 0,       0,      1]
]
where
    θ ∈ ℝ

# Rotation about Y-axis
mat_rot_y(θ) = [
    [cos(θ),  0, sin(θ), 0],
    [0,       1, 0,      0],
    [-sin(θ), 0, cos(θ), 0],
    [0,       0, 0,      1]
]
where
    θ ∈ ℝ

# Rotation about Z-axis
mat_rot_z(θ) = [
    [cos(θ),  sin(θ), 0, 0],
    [-sin(θ), cos(θ), 0, 0],
    [0,       0,      1, 0],
    [0,       0,      0, 1]
]
where
    θ ∈ ℝ

# Translation matrix
mat_trans(x, y, z) = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [x, y, z, 1]
]
where
    x ∈ ℝ
    y ∈ ℝ
    z ∈ ℝ

# Scale matrix
mat_scale(x, y, z) = [
    [x, 0, 0, 0],
    [0, y, 0, 0],
    [0, 0, z, 0],
    [0, 0, 0, 1]
]
where
    x ∈ ℝ
    y ∈ ℝ
    z ∈ ℝ

# Perspective projection (fov in degrees)
mat_proj(fov, aspect, near, far) =
    let f = 1 / tan(fov * 0.5 * π / 180) in
    let nf = far - near in
    [
        [aspect * f, 0, 0,                  0],
        [0,          f, 0,                  0],
        [0,          0, far / nf,           1],
        [0,          0, -far * near / nf,   0]
    ]
where
    fov ∈ ℝ
    aspect ∈ ℝ
    near ∈ ℝ
    far ∈ ℝ

# ═══════════════════════════════════════════════════════════════
# Matrix-Vector Multiplication: v' = v · M
# ═══════════════════════════════════════════════════════════════

# v'[j] = Σᵢ v[i] * M[i][j]
mat_mul_vec(M, v) = [
    ∑(i∈0..4) v[i] * M[i][0],
    ∑(i∈0..4) v[i] * M[i][1],
    ∑(i∈0..4) v[i] * M[i][2],
    ∑(i∈0..4) v[i] * M[i][3]
]
where
    M ∈ ℝ⁴ˣ⁴
    v ∈ ℝ⁴

# ═══════════════════════════════════════════════════════════════
# Matrix-Matrix Multiplication: C = A · B
# ═══════════════════════════════════════════════════════════════

# C[i][j] = Σₖ A[i][k] * B[k][j]
mat_mul(A, B) = [
    [∑(k∈0..4) A[0][k] * B[k][0], ∑(k∈0..4) A[0][k] * B[k][1],
     ∑(k∈0..4) A[0][k] * B[k][2], ∑(k∈0..4) A[0][k] * B[k][3]],
    [∑(k∈0..4) A[1][k] * B[k][0], ∑(k∈0..4) A[1][k] * B[k][1],
     ∑(k∈0..4) A[1][k] * B[k][2], ∑(k∈0..4) A[1][k] * B[k][3]],
    [∑(k∈0..4) A[2][k] * B[k][0], ∑(k∈0..4) A[2][k] * B[k][1],
     ∑(k∈0..4) A[2][k] * B[k][2], ∑(k∈0..4) A[2][k] * B[k][3]],
    [∑(k∈0..4) A[3][k] * B[k][0], ∑(k∈0..4) A[3][k] * B[k][1],
     ∑(k∈0..4) A[3][k] * B[k][2], ∑(k∈0..4) A[3][k] * B[k][3]]
]
where
    A ∈ ℝ⁴ˣ⁴
    B ∈ ℝ⁴ˣ⁴

# ═══════════════════════════════════════════════════════════════
# Camera Matrices
# ═══════════════════════════════════════════════════════════════

# Quick inverse for orthonormal rotation + translation
mat_qinv(M) =
    let m0 = mat_row3(M, 0) in
    let m1 = mat_row3(M, 1) in
    let m2 = mat_row3(M, 2) in
    let t = mat_row3(M, 3) in
    mat_from_rows(
        [m0.x, m1.x, m2.x, 0],
        [m0.y, m1.y, m2.y, 0],
        [m0.z, m1.z, m2.z, 0],
        [-vec_dot(t, m0), -vec_dot(t, m1), -vec_dot(t, m2), 1]
    )
where
    M ∈ ℝ⁴ˣ⁴

# Look-at: orthonormal basis from position/target/up
mat_point_at(pos, target, up) =
    let fwd = vec_sub(target, pos) in
    let fwd_len_sq = vec_dot(fwd, fwd) in
    if fwd_len_sq < ε then mat_ident() else
    let fwd_n = vec_norm(fwd) in
    let up_proj = vec_sub(up, vec_mul(fwd_n, vec_dot(up, fwd_n))) in
    let up_len_sq = vec_dot(up_proj, up_proj) in
    if up_len_sq < ε then mat_ident() else
    let up_n = vec_norm(up_proj) in
    let right = vec_cross(up_n, fwd_n) in
    mat_from_rows(
        [right.x, right.y, right.z, 0],
        [up_n.x, up_n.y, up_n.z, 0],
        [fwd_n.x, fwd_n.y, fwd_n.z, 0],
        [pos.x, pos.y, pos.z, 1]
    )
where
    pos ∈ ℝ⁴
    target ∈ ℝ⁴
    up ∈ ℝ⁴

# ═══════════════════════════════════════════════════════════════
# Geometry Operations
# ═══════════════════════════════════════════════════════════════

# Plane-line intersection
intersect_plane(n, d, start, end) =
    let ad = vec_dot(start, n) in
    let bd = vec_dot(end, n) in
    let denom = bd - ad in
    if |denom| < ε then start else
    let t = clamp((d - ad) / denom, 0, 1) in
    vec_lerp(start, end, t)
where
    n ∈ ℝ⁴
    d ∈ ℝ
    start ∈ ℝ⁴
    end ∈ ℝ⁴
