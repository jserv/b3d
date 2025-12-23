#!/usr/bin/env python3
"""
B3D Math DSL Code Generator - I❤LA Style

Generates C and LaTeX from mathematical notation.
Supports both Unicode symbols (∑, ∈, θ) and LaTeX escapes (\\sum, \\in, \\theta).

Usage:
    python3 scripts/gen-math.py --dsl src/math.dsl -o src/math-gen.inc
    python3 scripts/gen-math.py --dsl src/math.dsl --latex docs/math.tex
"""

import re
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional


class Mode(Enum):
    FLOAT = "float"
    FIXED = "fixed"


class TokenType(Enum):
    # Literals
    IDENT = auto()
    NUMBER = auto()
    STRING = auto()

    # Unicode math operators
    DOT = auto()  # · (middle dot)
    CROSS = auto()  # ×
    NORM_OPEN = auto()  # ‖
    NORM_CLOSE = auto()  # ‖
    SQRT = auto()  # √
    SUM = auto()  # ∑
    PROD = auto()  # ∏
    IN = auto()  # ∈
    TRANSPOSE = auto()  # ᵀ
    INVERSE = auto()  # ⁻¹

    # Greek letters
    THETA = auto()  # θ
    PI = auto()  # π
    EPSILON = auto()  # ε
    DELTA = auto()  # δ

    # Subscripts
    SUBSCRIPT = auto()  # ᵢⱼₖ₀₁₂₃

    # Type symbols
    REAL = auto()  # ℝ
    INT = auto()  # ℤ
    SUPERSCRIPT = auto()  # ³ ⁴ etc

    # Standard operators
    PLUS = auto()
    MINUS = auto()
    STAR = auto()
    SLASH = auto()
    CARET = auto()
    EQ = auto()
    LT = auto()
    GT = auto()
    LE = auto()
    GE = auto()
    PIPE = auto()  # |

    # Delimiters
    LPAREN = auto()
    RPAREN = auto()
    LBRACKET = auto()
    RBRACKET = auto()
    LBRACE = auto()
    RBRACE = auto()
    COMMA = auto()
    COLON = auto()
    SEMICOLON = auto()
    NEWLINE = auto()

    # Keywords
    WHERE = auto()
    LET = auto()
    IN_KW = auto()
    IF = auto()
    THEN = auto()
    ELSE = auto()
    SIN = auto()
    COS = auto()
    TAN = auto()
    ABS = auto()
    FLOOR = auto()
    MIN = auto()
    MAX = auto()
    CLAMP = auto()

    EOF = auto()


@dataclass
class Token:
    type: TokenType
    value: str
    line: int
    col: int


# Unicode character mappings
UNICODE_MAP = {
    "·": TokenType.DOT,
    "⋅": TokenType.DOT,
    "×": TokenType.CROSS,
    "‖": TokenType.NORM_OPEN,  # Context determines open/close
    "√": TokenType.SQRT,
    "∑": TokenType.SUM,
    "∏": TokenType.PROD,
    "∈": TokenType.IN,
    "ᵀ": TokenType.TRANSPOSE,
    "θ": TokenType.THETA,
    "π": TokenType.PI,
    "ε": TokenType.EPSILON,
    "δ": TokenType.DELTA,
    "ℝ": TokenType.REAL,
    "ℤ": TokenType.INT,
}

# LaTeX-like escape sequences (non-Unicode alternatives)
LATEX_ESCAPES = {
    # Operators
    "sum": TokenType.SUM,
    "prod": TokenType.PROD,
    "in": TokenType.IN,
    "cdot": TokenType.DOT,
    "dot": TokenType.DOT,
    "times": TokenType.CROSS,
    "cross": TokenType.CROSS,
    "sqrt": TokenType.SQRT,
    "norm": TokenType.NORM_OPEN,
    "T": TokenType.TRANSPOSE,
    "transpose": TokenType.TRANSPOSE,
    # Greek letters
    "theta": TokenType.THETA,
    "pi": TokenType.PI,
    "epsilon": TokenType.EPSILON,
    "eps": TokenType.EPSILON,
    "delta": TokenType.DELTA,
    # Type symbols
    "R": TokenType.REAL,
    "Real": TokenType.REAL,
    "Z": TokenType.INT,
    "Int": TokenType.INT,
}

SUBSCRIPTS = {
    "ᵢ": "i",
    "ⱼ": "j",
    "ₖ": "k",
    "ₗ": "l",
    "ₘ": "m",
    "ₙ": "n",
    "₀": "0",
    "₁": "1",
    "₂": "2",
    "₃": "3",
    "₄": "4",
}

SUPERSCRIPTS = {
    "⁰": "0",
    "¹": "1",
    "²": "2",
    "³": "3",
    "⁴": "4",
    "⁵": "5",
    "ˣ": "x",
}

KEYWORDS = {
    "where": TokenType.WHERE,
    "let": TokenType.LET,
    "in": TokenType.IN_KW,
    "if": TokenType.IF,
    "then": TokenType.THEN,
    "else": TokenType.ELSE,
    "sin": TokenType.SIN,
    "cos": TokenType.COS,
    "tan": TokenType.TAN,
    "abs": TokenType.ABS,
    "floor": TokenType.FLOOR,
    "min": TokenType.MIN,
    "max": TokenType.MAX,
    "clamp": TokenType.CLAMP,
}


class Lexer:
    """Unicode-aware lexer for I❤LA-style math DSL."""

    def __init__(self, source: str):
        self.source = source
        self.pos = 0
        self.line = 1
        self.col = 1
        self.tokens: list[Token] = []
        self.bracket_depth = 0  # Track nesting for multi-line expressions
        self.paren_depth = 0  # Track parentheses for multi-line expressions

    def peek(self, offset: int = 0) -> str:
        pos = self.pos + offset
        if pos >= len(self.source):
            return "\0"
        return self.source[pos]

    def advance(self) -> str:
        ch = self.peek()
        self.pos += 1
        if ch == "\n":
            self.line += 1
            self.col = 1
        else:
            self.col += 1
        return ch

    def skip_whitespace_and_comments(self):
        while self.pos < len(self.source):
            ch = self.peek()
            if ch == "#":
                while self.peek() not in ("\n", "\0"):
                    self.advance()
            elif ch in " \t\r":
                self.advance()
            else:
                break

    def read_number(self) -> Token:
        start_col = self.col
        result = []
        has_decimal = False
        while True:
            ch = self.peek()
            if ch.isdigit():
                result.append(self.advance())
            elif ch == "." and not has_decimal:
                # Don't consume if followed by another '.' (range operator ..)
                if self.peek(1) == ".":
                    break
                result.append(self.advance())
                has_decimal = True
            else:
                break
        return Token(TokenType.NUMBER, "".join(result), self.line, start_col)

    def read_ident(self) -> Token:
        start_col = self.col
        result = []

        # First character
        ch = self.advance()
        result.append(ch)

        # Rest of identifier
        while True:
            ch = self.peek()
            if ch.isalnum() or ch == "_":
                result.append(self.advance())
            elif ch in SUBSCRIPTS:
                # Convert subscript to bracket notation internally
                result.append("_")
                result.append(SUBSCRIPTS[self.advance()])
            else:
                break

        name = "".join(result)

        # Check for keywords
        if name.lower() in KEYWORDS:
            return Token(KEYWORDS[name.lower()], name, self.line, start_col)

        return Token(TokenType.IDENT, name, self.line, start_col)

    def tokenize(self) -> list[Token]:
        while self.pos < len(self.source):
            self.skip_whitespace_and_comments()
            if self.pos >= len(self.source):
                break

            ch = self.peek()
            start_line, start_col = self.line, self.col

            # Newlines (significant for where blocks, but ignored inside brackets/parens)
            if ch == "\n":
                self.advance()
                if self.bracket_depth == 0 and self.paren_depth == 0:
                    self.tokens.append(
                        Token(TokenType.NEWLINE, "\n", start_line, start_col)
                    )
                continue

            # LaTeX-like escape sequences: \sum, \in, \theta, etc.
            if ch == "\\":
                self.advance()
                # Read the escape name
                escape_name = []
                while self.peek().isalnum():
                    escape_name.append(self.advance())
                name = "".join(escape_name)

                if name in LATEX_ESCAPES:
                    self.tokens.append(
                        Token(LATEX_ESCAPES[name], "\\" + name, start_line, start_col)
                    )
                else:
                    # Unknown escape - treat as identifier
                    self.tokens.append(
                        Token(TokenType.IDENT, name, start_line, start_col)
                    )
                continue

            # Unicode operators
            if ch in UNICODE_MAP:
                self.advance()
                self.tokens.append(Token(UNICODE_MAP[ch], ch, start_line, start_col))
                continue

            # Subscripts as standalone tokens
            if ch in SUBSCRIPTS:
                self.advance()
                self.tokens.append(
                    Token(TokenType.SUBSCRIPT, SUBSCRIPTS[ch], start_line, start_col)
                )
                continue

            # Superscripts
            if ch in SUPERSCRIPTS:
                result = []
                while self.peek() in SUPERSCRIPTS:
                    result.append(SUPERSCRIPTS[self.advance()])
                self.tokens.append(
                    Token(TokenType.SUPERSCRIPT, "".join(result), start_line, start_col)
                )
                continue

            # Numbers
            if ch.isdigit():
                self.tokens.append(self.read_number())
                continue

            # Identifiers (including Greek letters as part of names)
            if ch.isalpha() or ch == "_" or ch in "αβγδεζηθικλμνξοπρστυφχψω":
                self.tokens.append(self.read_ident())
                continue

            # Two-character operators
            if ch == "." and self.peek(1) == ".":
                self.advance()
                self.advance()
                self.tokens.append(Token(TokenType.IDENT, "..", start_line, start_col))
                continue

            # Single period for field access
            if ch == ".":
                self.advance()
                self.tokens.append(Token(TokenType.DOT, ".", start_line, start_col))
                continue

            if ch == "-" and self.peek(1) == ">":
                self.advance()
                self.advance()
                self.tokens.append(Token(TokenType.IDENT, "->", start_line, start_col))
                continue

            # || as alternative to ‖ for norm brackets
            if ch == "|" and self.peek(1) == "|":
                self.advance()
                self.advance()
                self.tokens.append(
                    Token(TokenType.NORM_OPEN, "||", start_line, start_col)
                )
                continue

            # ^N for superscripts (e.g., ^4 instead of ⁴)
            if ch == "^" and self.peek(1).isdigit():
                self.advance()  # consume ^
                result = []
                while self.peek().isdigit() or self.peek() == "x":
                    result.append(self.advance())
                self.tokens.append(
                    Token(TokenType.SUPERSCRIPT, "".join(result), start_line, start_col)
                )
                continue

            # Single-character operators
            simple = {
                "+": TokenType.PLUS,
                "-": TokenType.MINUS,
                "*": TokenType.STAR,
                "/": TokenType.SLASH,
                "^": TokenType.CARET,
                "=": TokenType.EQ,
                "<": TokenType.LT,
                ">": TokenType.GT,
                "|": TokenType.PIPE,
                "(": TokenType.LPAREN,
                ")": TokenType.RPAREN,
                "[": TokenType.LBRACKET,
                "]": TokenType.RBRACKET,
                "{": TokenType.LBRACE,
                "}": TokenType.RBRACE,
                ",": TokenType.COMMA,
                ":": TokenType.COLON,
                ";": TokenType.SEMICOLON,
            }

            if ch in simple:
                self.advance()
                # Track nesting depth for multi-line expression support
                if ch == "[":
                    self.bracket_depth += 1
                elif ch == "]":
                    self.bracket_depth = max(0, self.bracket_depth - 1)
                elif ch == "(":
                    self.paren_depth += 1
                elif ch == ")":
                    self.paren_depth = max(0, self.paren_depth - 1)
                self.tokens.append(Token(simple[ch], ch, start_line, start_col))
                continue

            # Skip unknown characters
            self.advance()

        self.tokens.append(Token(TokenType.EOF, "", self.line, self.col))
        return self.tokens


# AST Nodes
@dataclass
class Expr:
    pass


@dataclass
class NumExpr(Expr):
    value: str


@dataclass
class VarExpr(Expr):
    name: str


@dataclass
class BinOpExpr(Expr):
    op: str
    left: Expr
    right: Expr


@dataclass
class UnaryExpr(Expr):
    op: str
    operand: Expr


@dataclass
class CallExpr(Expr):
    func: str
    args: list[Expr]


@dataclass
class IndexExpr(Expr):
    base: Expr
    index: Expr


@dataclass
class MatrixIndexExpr(Expr):
    base: Expr
    row: Expr
    col: Expr


@dataclass
class DotAccessExpr(Expr):
    base: Expr
    field: str


@dataclass
class SumExpr(Expr):
    var: str
    range_expr: str  # "xyz" or "0..4"
    body: Expr


@dataclass
class NormExpr(Expr):
    operand: Expr


@dataclass
class VectorExpr(Expr):
    elements: list[Expr]


@dataclass
class ComprehensionExpr(Expr):
    body: Expr
    var: str
    range_expr: str


@dataclass
class MatrixExpr(Expr):
    rows: list[list[Expr]]


@dataclass
class LetExpr(Expr):
    bindings: list[tuple[str, Expr]]
    body: Expr


@dataclass
class IfExpr(Expr):
    cond: Expr
    then_expr: Expr
    else_expr: Expr


@dataclass
class Param:
    name: str
    type_str: str  # "ℝ⁴" or "scalar" etc


@dataclass
class FuncDef:
    name: str
    params: list[Param]
    return_type: str
    body: Expr


class Parser:
    """Recursive descent parser for I❤LA-style DSL."""

    def __init__(self, tokens: list[Token]):
        self.tokens = tokens
        self.pos = 0

    def peek(self, offset: int = 0) -> Token:
        pos = self.pos + offset
        if pos >= len(self.tokens):
            return self.tokens[-1]  # EOF
        return self.tokens[pos]

    def advance(self) -> Token:
        tok = self.peek()
        self.pos += 1
        return tok

    def expect(self, token_type: TokenType) -> Token:
        tok = self.advance()
        if tok.type != token_type:
            raise SyntaxError(
                f"Expected {token_type}, got {tok.type} ({tok.value!r}) "
                f"at line {tok.line}, col {tok.col}"
            )
        return tok

    def skip_newlines(self):
        while self.peek().type == TokenType.NEWLINE:
            self.advance()

    def parse(self) -> list[FuncDef]:
        funcs = []
        while self.peek().type != TokenType.EOF:
            self.skip_newlines()
            if self.peek().type == TokenType.EOF:
                break
            func = self.parse_func()
            if func:
                funcs.append(func)
        return funcs

    def parse_func(self) -> Optional[FuncDef]:
        """Parse: name(params) = expr where types"""
        self.skip_newlines()

        # Function name
        if self.peek().type != TokenType.IDENT:
            self.advance()  # Skip unknown
            return None

        name = self.advance().value

        # Parameters
        self.expect(TokenType.LPAREN)
        params = self.parse_param_list()
        self.expect(TokenType.RPAREN)

        # Equals sign
        self.expect(TokenType.EQ)

        # Body expression
        body = self.parse_expr()

        # Optional where clause
        self.skip_newlines()
        param_types = {}
        if self.peek().type == TokenType.WHERE:
            self.advance()
            self.skip_newlines()
            param_types = self.parse_where_block()

        # Build params with types
        typed_params = []
        for p in params:
            type_str = param_types.get(p, "scalar")
            typed_params.append(Param(p, type_str))

        # Infer return type from body
        return_type = self.infer_return_type(body, param_types)

        return FuncDef(name, typed_params, return_type, body)

    def parse_param_list(self) -> list[str]:
        params = []
        if self.peek().type == TokenType.RPAREN:
            return params

        # Handle both IDENT and Greek letters (THETA, etc.)
        tok = self.advance()
        if tok.type == TokenType.THETA:
            params.append("a")  # θ → a (angle) in C
        else:
            params.append(tok.value)

        while self.peek().type == TokenType.COMMA:
            self.advance()
            tok = self.advance()
            if tok.type == TokenType.THETA:
                params.append("a")
            else:
                params.append(tok.value)
        return params

    def parse_where_block(self) -> dict[str, str]:
        """Parse type declarations in where block."""
        types = {}
        while True:
            self.skip_newlines()
            tok = self.peek()

            # Check if this looks like a new function definition
            if tok.type == TokenType.IDENT and self.peek(1).type == TokenType.LPAREN:
                break

            # Handle both IDENT and Greek letters
            if tok.type == TokenType.IDENT:
                name = self.advance().value
            elif tok.type == TokenType.THETA:
                self.advance()
                name = "a"  # θ → a
            else:
                break

            if self.peek().type != TokenType.IN:
                self.pos -= 1
                break
            self.advance()  # consume ∈

            type_str = self.parse_type()
            types[name] = type_str

        return types

    def parse_type(self) -> str:
        """Parse type: ℝ, ℝ³, ℝ⁴, ℝ⁴ˣ⁴, ℤ"""
        result = []

        if self.peek().type == TokenType.REAL:
            self.advance()
            result.append("ℝ")
        elif self.peek().type == TokenType.INT:
            self.advance()
            result.append("ℤ")
        else:
            return "scalar"

        # Superscripts for dimensions
        if self.peek().type == TokenType.SUPERSCRIPT:
            result.append(self.advance().value)

        return "".join(result)

    def infer_return_type(self, body: Expr, types: dict) -> str:
        """Infer return type from expression structure."""
        if isinstance(body, MatrixExpr):
            return "mat4"
        if isinstance(body, VectorExpr) or isinstance(body, ComprehensionExpr):
            return "vec4"
        if isinstance(body, SumExpr):
            return "scalar"
        if isinstance(body, NormExpr):
            return "scalar"
        if isinstance(body, IfExpr):
            # Check both branches, prefer then_expr but try else_expr too
            then_type = self.infer_return_type(body.then_expr, types)
            if then_type != "scalar":
                return then_type
            return self.infer_return_type(body.else_expr, types)
        if isinstance(body, LetExpr):
            return self.infer_return_type(body.body, types)
        if isinstance(body, VarExpr):
            # Look up variable type from types dict
            var_type = types.get(body.name, "scalar")
            if var_type in ("ℝ⁴", "ℝ³"):
                return "vec4"
            if var_type in ("ℝ⁴ˣ⁴",):
                return "mat4"
        if isinstance(body, CallExpr):
            if body.func in ("vec_dot", "vec_length", "vec_length_sq"):
                return "scalar"
            if body.func.startswith("mat_"):
                if "vec" in body.func:
                    return "vec4"
                return "mat4"
            if body.func.startswith("vec_"):
                return "vec4"
        return "scalar"

    def parse_expr(self) -> Expr:
        self.skip_newlines()
        return self.parse_let()

    def parse_let(self) -> Expr:
        """Parse let bindings."""
        self.skip_newlines()
        if self.peek().type == TokenType.LET:
            self.advance()
            bindings = []
            name = self.advance().value
            self.expect(TokenType.EQ)
            value = self.parse_conditional()
            bindings.append((name, value))

            while self.peek().type == TokenType.SEMICOLON:
                self.advance()
                self.skip_newlines()
                if self.peek().type == TokenType.LET:
                    self.advance()
                name = self.advance().value
                self.expect(TokenType.EQ)
                value = self.parse_conditional()
                bindings.append((name, value))

            if self.peek().type == TokenType.IN_KW:
                self.advance()
                body = self.parse_expr()
            else:
                body = self.parse_expr()

            return LetExpr(bindings, body)

        return self.parse_conditional()

    def parse_conditional(self) -> Expr:
        """Parse if/then/else."""
        if self.peek().type == TokenType.IF:
            self.advance()
            cond = self.parse_comparison()
            self.skip_newlines()
            self.expect(TokenType.THEN)
            # Both branches can contain let expressions
            then_expr = self.parse_let()
            self.skip_newlines()
            self.expect(TokenType.ELSE)
            else_expr = self.parse_let()
            return IfExpr(cond, then_expr, else_expr)

        return self.parse_comparison()

    def parse_comparison(self) -> Expr:
        left = self.parse_additive()

        while self.peek().type in (
            TokenType.LT,
            TokenType.GT,
            TokenType.LE,
            TokenType.GE,
        ):
            op = self.advance().value
            right = self.parse_additive()
            left = BinOpExpr(op, left, right)

        return left

    def parse_additive(self) -> Expr:
        left = self.parse_multiplicative()

        while self.peek().type in (TokenType.PLUS, TokenType.MINUS):
            op = self.advance().value
            right = self.parse_multiplicative()
            left = BinOpExpr(op, left, right)

        return left

    def parse_multiplicative(self) -> Expr:
        left = self.parse_unary()

        while self.peek().type in (
            TokenType.STAR,
            TokenType.SLASH,
            TokenType.DOT,
            TokenType.CROSS,
        ):
            # Don't consume DOT if followed by field access (ident)
            if self.peek().type == TokenType.DOT:
                # Check if this is field access (dot followed by ident)
                if self.peek(1).type == TokenType.IDENT:
                    break
            op_tok = self.advance()
            op = op_tok.value
            if op_tok.type == TokenType.DOT:
                op = "dot"
            elif op_tok.type == TokenType.CROSS:
                op = "cross"
            right = self.parse_unary()
            left = BinOpExpr(op, left, right)

        return left

    def parse_unary(self) -> Expr:
        if self.peek().type == TokenType.MINUS:
            self.advance()
            operand = self.parse_unary()
            return UnaryExpr("-", operand)

        if self.peek().type == TokenType.SQRT:
            self.advance()
            if self.peek().type == TokenType.LPAREN:
                self.advance()
                operand = self.parse_expr()
                self.expect(TokenType.RPAREN)
            else:
                operand = self.parse_primary()
            return CallExpr("sqrt", [operand])

        return self.parse_postfix()

    def parse_postfix(self) -> Expr:
        left = self.parse_primary()

        while True:
            if self.peek().type == TokenType.LBRACKET:
                self.advance()
                index = self.parse_expr()
                self.expect(TokenType.RBRACKET)

                # Check for second index
                if self.peek().type == TokenType.LBRACKET:
                    self.advance()
                    index2 = self.parse_expr()
                    self.expect(TokenType.RBRACKET)
                    left = MatrixIndexExpr(left, index, index2)
                else:
                    left = IndexExpr(left, index)

            elif (
                self.peek().type == TokenType.DOT
                and self.peek(1).type == TokenType.IDENT
            ):
                self.advance()
                field = self.advance().value
                left = DotAccessExpr(left, field)

            elif self.peek().type == TokenType.TRANSPOSE:
                self.advance()
                left = UnaryExpr("T", left)

            else:
                break

        return left

    def parse_primary(self) -> Expr:
        tok = self.peek()

        # Number
        if tok.type == TokenType.NUMBER:
            self.advance()
            return NumExpr(tok.value)

        # Sum: ∑(i∈range) body
        if tok.type == TokenType.SUM:
            self.advance()
            self.expect(TokenType.LPAREN)
            var = self.advance().value
            self.expect(TokenType.IN)
            range_expr = self.parse_range()
            self.expect(TokenType.RPAREN)
            body = self.parse_multiplicative()
            return SumExpr(var, range_expr, body)

        # Norm: ‖expr‖
        if tok.type == TokenType.NORM_OPEN:
            self.advance()
            operand = self.parse_expr()
            self.expect(TokenType.NORM_OPEN)  # Same token for close
            return NormExpr(operand)

        # Absolute value: |expr|
        if tok.type == TokenType.PIPE:
            self.advance()
            operand = self.parse_expr()
            self.expect(TokenType.PIPE)
            return CallExpr("abs", [operand])

        # Vector literal or comprehension: [...]
        if tok.type == TokenType.LBRACKET:
            return self.parse_vector_or_matrix()

        # Parenthesized expression
        if tok.type == TokenType.LPAREN:
            self.advance()
            expr = self.parse_expr()
            self.expect(TokenType.RPAREN)
            return expr

        # Built-in functions
        if tok.type in (
            TokenType.SIN,
            TokenType.COS,
            TokenType.TAN,
            TokenType.ABS,
            TokenType.FLOOR,
            TokenType.MIN,
            TokenType.MAX,
            TokenType.CLAMP,
        ):
            func = self.advance().value.lower()
            if self.peek().type == TokenType.LPAREN:
                self.advance()
                args = self.parse_arg_list()
                self.expect(TokenType.RPAREN)
            else:
                args = [self.parse_unary()]
            return CallExpr(func, args)

        # Constants
        if tok.type == TokenType.PI:
            self.advance()
            return VarExpr("PI")

        if tok.type == TokenType.EPSILON:
            self.advance()
            return VarExpr("EPSILON")

        if tok.type == TokenType.THETA:
            self.advance()
            return VarExpr("a")  # θ → a (angle) in C

        if tok.type == TokenType.DELTA:
            self.advance()
            # Kronecker delta - need subscripts
            if self.peek().type == TokenType.SUBSCRIPT:
                i = self.advance().value
                if self.peek().type == TokenType.SUBSCRIPT:
                    j = self.advance().value
                    return CallExpr("kronecker", [VarExpr(i), VarExpr(j)])
            return VarExpr("delta")

        # Identifier (variable or function call)
        if tok.type == TokenType.IDENT:
            name = self.advance().value

            # Function call - check before subscript processing
            if self.peek().type == TokenType.LPAREN:
                self.advance()
                args = self.parse_arg_list()
                self.expect(TokenType.RPAREN)
                return CallExpr(name, args)

            # Check for subscript notation: v_i (single letter base only)
            # Only applies to patterns like a_i, M_0, not vec_dot or up_n
            # Subscript suffixes are limited to typical indices: i, j, k, l, m, n, 0-9
            # But NOT when they look like variable suffixes (up_n means normalized up)
            if "_" in name:
                parts = name.split("_")
                # Only treat as subscript if:
                # 1. Base is single char (a, v, M) - NOT 2 chars like "up"
                # 2. All suffixes are single digit or typical loop index
                valid_indices = {"i", "j", "k", "l", "m", "0", "1", "2", "3", "4"}
                if len(parts[0]) == 1 and all(p in valid_indices for p in parts[1:]):
                    base = VarExpr(parts[0])
                    for idx in parts[1:]:
                        base = IndexExpr(base, VarExpr(idx))
                    return base

            return VarExpr(name)

        # Unknown - skip
        self.advance()
        return VarExpr("_unknown_")

    def parse_range(self) -> str:
        """Parse range like xyz or 0..4"""
        parts = []
        tok = self.advance()
        parts.append(tok.value)

        if self.peek().type == TokenType.IDENT and self.peek().value == "..":
            self.advance()  # ..
            parts.append("..")
            parts.append(self.advance().value)

        return "".join(parts)

    def parse_arg_list(self) -> list[Expr]:
        args = []
        if self.peek().type == TokenType.RPAREN:
            return args

        args.append(self.parse_expr())
        while self.peek().type == TokenType.COMMA:
            self.advance()
            args.append(self.parse_expr())
        return args

    def parse_vector_or_matrix(self) -> Expr:
        """Parse [...] - vector, comprehension, or matrix."""
        self.expect(TokenType.LBRACKET)

        # Check for matrix: [[...], [...]]
        if self.peek().type == TokenType.LBRACKET:
            rows = []
            while self.peek().type == TokenType.LBRACKET:
                self.advance()
                row = []
                row.append(self.parse_expr())
                while self.peek().type == TokenType.COMMA:
                    self.advance()
                    if self.peek().type == TokenType.RBRACKET:
                        break
                    row.append(self.parse_expr())
                self.expect(TokenType.RBRACKET)
                rows.append(row)
                if self.peek().type == TokenType.COMMA:
                    self.advance()
                if self.peek().type == TokenType.RBRACKET:
                    break
            self.expect(TokenType.RBRACKET)
            return MatrixExpr(rows)

        # First element
        first = self.parse_expr()

        # Comprehension: [expr | var∈range]
        if self.peek().type == TokenType.PIPE:
            self.advance()
            var = self.advance().value
            self.expect(TokenType.IN)
            range_expr = self.parse_range()
            self.expect(TokenType.RBRACKET)
            return ComprehensionExpr(first, var, range_expr)

        # Vector literal: [a, b, c, d] or could be matrix rows
        elements = [first]
        while self.peek().type == TokenType.COMMA:
            self.advance()
            if self.peek().type == TokenType.RBRACKET:
                break
            elements.append(self.parse_expr())

        self.expect(TokenType.RBRACKET)

        # If first element is a vector, this is a matrix
        if isinstance(first, VectorExpr):
            return MatrixExpr(
                [e.elements if isinstance(e, VectorExpr) else [e] for e in elements]
            )

        return VectorExpr(elements)


# ==============================================================================
# Code Generator
# ==============================================================================


class CodeGen:
    """Generate C code from AST."""

    def __init__(self, mode: Mode, suffix: str = ""):
        self.mode = mode
        self.suffix = suffix
        self.emit_int_context = False  # Whether to emit integers without .0f

    def c_type(self, type_str: str) -> str:
        if type_str in ("ℝ", "scalar", "ℝ¹"):
            return "float"
        if type_str in ("ℝ³", "ℝ⁴", "ℝ4", "vec4", "vec3"):
            return "b3d_vec_t"
        if type_str in ("ℝ⁴ˣ⁴", "ℝ4x4", "mat4"):
            return "b3d_mat_t"
        if type_str in ("ℤ", "int"):
            return "int"
        # Default based on name pattern
        return "float"

    def has_nested_conditional_lets(self, expr: Expr) -> bool:
        """Check if expression has let bindings inside conditionals."""
        if isinstance(expr, LetExpr):
            return self.has_nested_conditional_lets(expr.body)
        if isinstance(expr, IfExpr):
            # Check if either branch contains let bindings
            then_has = isinstance(
                expr.then_expr, LetExpr
            ) or self.has_nested_conditional_lets(expr.then_expr)
            else_has = isinstance(
                expr.else_expr, LetExpr
            ) or self.has_nested_conditional_lets(expr.else_expr)
            return then_has or else_has
        return False

    def emit_complex_body(self, expr: Expr, ret_type: str, indent: int = 1) -> str:
        """Emit complex function body with proper if/else blocks for nested lets."""
        lines = []
        ind = "    " * indent

        # First, extract top-level let bindings
        while isinstance(expr, LetExpr):
            for name, value in expr.bindings:
                val_code = self.emit_expr(value)
                val_type = self.infer_c_type(value)
                lines.append(f"{ind}{val_type} {name} = {val_code};")
            expr = expr.body

        # Handle conditional with nested lets
        if isinstance(expr, IfExpr):
            cond = self.emit_expr(expr.cond)

            # Check if branches have let bindings
            then_has_lets = isinstance(
                expr.then_expr, LetExpr
            ) or self.has_nested_conditional_lets(expr.then_expr)
            else_has_lets = isinstance(
                expr.else_expr, LetExpr
            ) or self.has_nested_conditional_lets(expr.else_expr)

            if then_has_lets or else_has_lets:
                # Emit as if/else block
                lines.append(f"{ind}if ({cond}) {{")
                then_body = self.emit_complex_body(expr.then_expr, ret_type, indent + 1)
                lines.append(then_body)
                lines.append(f"{ind}}} else {{")
                else_body = self.emit_complex_body(expr.else_expr, ret_type, indent + 1)
                lines.append(else_body)
                lines.append(f"{ind}}}")
            else:
                # Simple ternary
                then_code = self.emit_expr(expr.then_expr)
                else_code = self.emit_expr(expr.else_expr)
                lines.append(f"{ind}return (({cond}) ? {then_code} : {else_code});")
        else:
            # Simple expression
            code = self.emit_expr(expr)
            lines.append(f"{ind}return {code};")

        return "\n".join(lines)

    def emit_func(self, func: FuncDef) -> str:
        ret_type = self.c_type(func.return_type)
        params = ", ".join(f"{self.c_type(p.type_str)} {p.name}" for p in func.params)
        func_name = f"b3d_{func.name}{self.suffix}"

        # Check if this function has complex nested conditionals with let bindings
        if self.has_nested_conditional_lets(func.body):
            body_lines = self.emit_complex_body(func.body, ret_type)
            return f"""static inline {ret_type} {func_name}({params})
{{
{body_lines}
}}"""

        # Simple case: extract let bindings
        declarations, final_expr = self.extract_lets(func.body)

        if declarations:
            decl_str = "\n".join(f"    {d}" for d in declarations)
            return f"""static inline {ret_type} {func_name}({params})
{{
{decl_str}

    return {final_expr};
}}"""
        else:
            body_code = self.emit_expr(func.body)
            return f"""static inline {ret_type} {func_name}({params})
{{
    return {body_code};
}}"""

    def extract_lets(self, expr: Expr) -> tuple[list[str], str]:
        """Extract let bindings into C variable declarations."""
        declarations = []

        while isinstance(expr, LetExpr):
            for name, value in expr.bindings:
                val_code = self.emit_expr(value)
                val_type = self.infer_c_type(value)
                declarations.append(f"{val_type} {name} = {val_code};")
            expr = expr.body

        final = self.emit_expr(expr)
        return declarations, final

    def infer_c_type(self, expr: Expr) -> str:
        if isinstance(expr, NormExpr):
            return "float"
        if isinstance(expr, CallExpr):
            if expr.func in ("vec_dot", "vec_length", "vec_length_sq", "sqrt"):
                return "float"
            # mat_row3 returns a vector, not a matrix
            if expr.func == "mat_row3":
                return "b3d_vec_t"
            if expr.func.startswith("vec_"):
                return "b3d_vec_t"
            if expr.func.startswith("mat_"):
                return "b3d_mat_t"
        if isinstance(expr, VectorExpr) or isinstance(expr, ComprehensionExpr):
            return "b3d_vec_t"
        if isinstance(expr, MatrixExpr):
            return "b3d_mat_t"
        return "float"

    def emit_expr(self, expr: Expr) -> str:
        if isinstance(expr, NumExpr):
            val = expr.value
            # Check if this should be emitted as integer
            if self.emit_int_context:
                # Strip any float suffix and return as int
                if "." in val:
                    val = val.split(".")[0]
                return val
            if "." not in val:
                val += ".0"
            return val + "f"

        if isinstance(expr, VarExpr):
            return self.emit_var(expr.name)

        if isinstance(expr, BinOpExpr):
            left = self.emit_expr(expr.left)
            right = self.emit_expr(expr.right)
            return self.emit_binop(expr.op, left, right)

        if isinstance(expr, UnaryExpr):
            operand = self.emit_expr(expr.operand)
            if expr.op == "-":
                return f"-({operand})"
            if expr.op == "T":
                return f"b3d_mat_transpose{self.suffix}({operand})"
            return f"{expr.op}({operand})"

        if isinstance(expr, CallExpr):
            args = self.emit_call_args(expr.func, expr.args)
            return self.emit_call(expr.func, args)

        if isinstance(expr, IndexExpr):
            base = self.emit_expr(expr.base)
            idx = self.emit_expr(expr.index)
            # Convert numeric index to component
            if idx in ("0", "0.0f"):
                return f"{base}.x"
            if idx in ("1", "1.0f"):
                return f"{base}.y"
            if idx in ("2", "2.0f"):
                return f"{base}.z"
            if idx in ("3", "3.0f"):
                return f"{base}.w"
            if idx in ("i", "j", "k"):
                return f"{base}[{idx}]"
            # Component access
            if idx in ("x", "y", "z", "w"):
                return f"{base}.{idx}"
            return f"{base}.{idx}"

        if isinstance(expr, MatrixIndexExpr):
            base = self.emit_expr(expr.base)
            row = self.emit_expr(expr.row)
            col = self.emit_expr(expr.col)
            # Remove .0f suffix from integer indices
            row = re.sub(r"\.0f$", "", row)
            col = re.sub(r"\.0f$", "", col)
            return f"{base}.m[{row}][{col}]"

        if isinstance(expr, DotAccessExpr):
            base = self.emit_expr(expr.base)
            return f"{base}.{expr.field}"

        if isinstance(expr, SumExpr):
            return self.emit_sum(expr)

        if isinstance(expr, NormExpr):
            operand = self.emit_expr(expr.operand)
            return f"b3d_vec_length{self.suffix}({operand})"

        if isinstance(expr, VectorExpr):
            elems = [self.emit_expr(e) for e in expr.elements]
            return "(b3d_vec_t){" + ", ".join(elems) + "}"

        if isinstance(expr, ComprehensionExpr):
            return self.emit_comprehension(expr)

        if isinstance(expr, MatrixExpr):
            rows = []
            for row in expr.rows:
                elems = [self.emit_expr(e) for e in row]
                rows.append("{" + ", ".join(elems) + "}")
            return "(b3d_mat_t){.m = {" + ", ".join(rows) + "}}"

        if isinstance(expr, LetExpr):
            # For inline use (shouldn't happen if extract_lets works)
            return self.emit_expr(expr.body)

        if isinstance(expr, IfExpr):
            cond = self.emit_expr(expr.cond)
            then = self.emit_expr(expr.then_expr)
            else_ = self.emit_expr(expr.else_expr)
            return f"(({cond}) ? {then} : {else_})"

        return "_unknown_"

    def emit_var(self, name: str) -> str:
        constants = {
            "PI": "B3D_PI" if self.mode == Mode.FLOAT else "B3D_FP_PI",
            "EPSILON": "B3D_EPSILON" if self.mode == Mode.FLOAT else "B3D_FP_EPSILON",
            "ZERO": "0.0f" if self.mode == Mode.FLOAT else "0",
            "ONE": "1.0f" if self.mode == Mode.FLOAT else "B3D_FP_ONE",
            "I": f"b3d_mat_ident{self.suffix}()",
        }
        return constants.get(name, name)

    def emit_binop(self, op: str, left: str, right: str) -> str:
        if self.mode == Mode.FLOAT:
            if op == "dot":
                return f"b3d_vec_dot{self.suffix}({left}, {right})"
            if op == "cross":
                return f"b3d_vec_cross{self.suffix}({left}, {right})"
            return f"(({left}) {op} ({right}))"
        else:
            # Fixed-point mode
            ops = {
                "+": f"B3D_FP_ADD({left}, {right})",
                "-": f"B3D_FP_SUB({left}, {right})",
                "*": f"B3D_FP_MUL({left}, {right})",
                "/": f"B3D_FP_DIV({left}, {right})",
                "dot": f"b3d_vec_dot{self.suffix}({left}, {right})",
                "cross": f"b3d_vec_cross{self.suffix}({left}, {right})",
            }
            return ops.get(op, f"(({left}) {op} ({right}))")

    # Map function names to parameter types that should be emitted as integers
    # Format: func_name -> list of indices that are integers
    INT_PARAM_FUNCS = {
        "mat_row3": [1],  # second param (row index) is int
    }

    def emit_call_args(self, func: str, args: list[Expr]) -> list[str]:
        """Emit function arguments, respecting integer parameter types."""
        int_indices = self.INT_PARAM_FUNCS.get(func, [])
        result = []
        for i, arg in enumerate(args):
            if i in int_indices:
                self.emit_int_context = True
                result.append(self.emit_expr(arg))
                self.emit_int_context = False
            else:
                result.append(self.emit_expr(arg))
        return result

    def emit_call(self, func: str, args: list[str]) -> str:
        # Math functions
        if self.mode == Mode.FLOAT:
            builtins = {
                "sqrt": lambda a: f"sqrtf({a[0]})",
                "sin": lambda a: f"sinf({a[0]})",
                "cos": lambda a: f"cosf({a[0]})",
                "tan": lambda a: f"tanf({a[0]})",
                "abs": lambda a: f"fabsf({a[0]})",
                "floor": lambda a: f"floorf({a[0]})",
                "min": lambda a: f"fminf({a[0]}, {a[1]})",
                "max": lambda a: f"fmaxf({a[0]}, {a[1]})",
                "clamp": lambda a: f"fminf(fmaxf({a[0]}, {a[1]}), {a[2]})",
                "kronecker": lambda a: f"(({a[0]}) == ({a[1]}) ? 1.0f : 0.0f)",
            }
        else:
            builtins = {
                "sqrt": lambda a: f"b3d_fp_sqrt({a[0]})",
                "sin": lambda a: f"b3d_fp_sin({a[0]})",
                "cos": lambda a: f"b3d_fp_cos({a[0]})",
                "tan": lambda a: f"b3d_fp_tan({a[0]})",
                "abs": lambda a: f"b3d_fp_abs({a[0]})",
                "floor": lambda a: f"B3D_FP_FLOOR({a[0]})",
                "min": lambda a: f"B3D_FP_MIN({a[0]}, {a[1]})",
                "max": lambda a: f"B3D_FP_MAX({a[0]}, {a[1]})",
                "clamp": lambda a: f"B3D_FP_MIN(B3D_FP_MAX({a[0]}, {a[1]}), {a[2]})",
                "kronecker": lambda a: f"(({a[0]}) == ({a[1]}) ? B3D_FP_ONE : 0)",
            }

        if func in builtins:
            return builtins[func](args)

        # B3D functions
        return f"b3d_{func}{self.suffix}({', '.join(args)})"

    def emit_sum(self, expr: SumExpr) -> str:
        """Expand sum to explicit additions."""
        indices = self.parse_range(expr.range_expr)
        terms = []

        for idx in indices:
            body_code = self.emit_expr(expr.body)
            # Substitute index variable
            body_code = self.substitute_index(body_code, expr.var, idx)
            terms.append(f"({body_code})")

        return " + ".join(terms)

    def emit_comprehension(self, expr: ComprehensionExpr) -> str:
        """Expand comprehension to vector literal."""
        indices = self.parse_range(expr.range_expr)
        elements = []

        for idx in indices:
            elem_code = self.emit_expr(expr.body)
            elem_code = self.substitute_index(elem_code, expr.var, idx)
            elements.append(elem_code)

        return "(b3d_vec_t){" + ", ".join(elements) + "}"

    def parse_range(self, range_str: str) -> list[str]:
        if range_str == "xyz":
            return ["x", "y", "z"]
        if range_str == "xyzw":
            return ["x", "y", "z", "w"]
        if ".." in range_str:
            start, end = range_str.split("..")
            return [str(i) for i in range(int(start), int(end))]
        return [range_str]

    def substitute_index(self, code: str, var: str, val: str) -> str:
        """Replace index variable with concrete value."""
        # Handle array access patterns
        code = re.sub(rf"\b{var}\b", val, code)

        # Fix component access for vector indices
        if val in ("x", "y", "z", "w"):
            # Replace [x] with .x for vectors (but not for matrices with second index)
            code = re.sub(r"\[" + val + r"\](?!\[)", f".{val}", code)
        elif val in ("0", "1", "2", "3"):
            component = ["x", "y", "z", "w"][int(val)]
            # Replace [0.0f] with [0] for matrix indices
            code = re.sub(r"\[" + val + r"\.0f\]", f"[{val}]", code)
            # Replace [0] with .x for vectors (only if not followed by another [)
            code = re.sub(r"(?<!\])\[" + val + r"\](?!\[)", f".{component}", code)

        # Clean up any remaining .0f in matrix indices
        code = re.sub(r"\.m\[(\d+)\.0f\]", r".m[\1]", code)

        return code


# ==============================================================================
# LaTeX Generator
# ==============================================================================


class LaTeXGen:
    """Generate LaTeX from AST."""

    def __init__(self):
        pass

    def emit_func(self, func: FuncDef) -> str:
        params = ", ".join(p.name for p in func.params)
        body = self.emit_expr(func.body)

        latex = f"\\textbf{{{func.name}}}({params}) &= {body}\n"
        return latex

    def emit_expr(self, expr: Expr) -> str:
        if isinstance(expr, NumExpr):
            return expr.value

        if isinstance(expr, VarExpr):
            return self.emit_var(expr.name)

        if isinstance(expr, BinOpExpr):
            left = self.emit_expr(expr.left)
            right = self.emit_expr(expr.right)
            return self.emit_binop(expr.op, left, right)

        if isinstance(expr, UnaryExpr):
            operand = self.emit_expr(expr.operand)
            if expr.op == "-":
                return f"-{operand}"
            if expr.op == "T":
                return f"{operand}^T"
            return f"{expr.op}({operand})"

        if isinstance(expr, CallExpr):
            args = [self.emit_expr(a) for a in expr.args]
            return self.emit_call(expr.func, args)

        if isinstance(expr, IndexExpr):
            base = self.emit_expr(expr.base)
            idx = self.emit_expr(expr.index)
            return f"{base}_{{{idx}}}"

        if isinstance(expr, MatrixIndexExpr):
            base = self.emit_expr(expr.base)
            row = self.emit_expr(expr.row)
            col = self.emit_expr(expr.col)
            return f"{base}_{{{row},{col}}}"

        if isinstance(expr, DotAccessExpr):
            base = self.emit_expr(expr.base)
            return f"{base}.{expr.field}"

        if isinstance(expr, SumExpr):
            body = self.emit_expr(expr.body)
            return f"\\sum_{{{expr.var} \\in {expr.range_expr}}} {body}"

        if isinstance(expr, NormExpr):
            operand = self.emit_expr(expr.operand)
            return f"\\|{operand}\\|"

        if isinstance(expr, VectorExpr):
            elems = [self.emit_expr(e) for e in expr.elements]
            return "\\begin{bmatrix}" + " \\\\ ".join(elems) + "\\end{bmatrix}"

        if isinstance(expr, ComprehensionExpr):
            body = self.emit_expr(expr.body)
            return f"[{body} \\mid {expr.var} \\in {expr.range_expr}]"

        if isinstance(expr, MatrixExpr):
            rows = []
            for row in expr.rows:
                elems = [self.emit_expr(e) for e in row]
                rows.append(" & ".join(elems))
            return "\\begin{bmatrix}" + " \\\\ ".join(rows) + "\\end{bmatrix}"

        if isinstance(expr, LetExpr):
            bindings = []
            for name, value in expr.bindings:
                bindings.append(f"{name} = {self.emit_expr(value)}")
            body = self.emit_expr(expr.body)
            return f"\\text{{let }} {', '.join(bindings)} \\text{{ in }} {body}"

        if isinstance(expr, IfExpr):
            cond = self.emit_expr(expr.cond)
            then = self.emit_expr(expr.then_expr)
            else_ = self.emit_expr(expr.else_expr)
            return f"\\begin{{cases}} {then} & \\text{{if }} {cond} \\\\ {else_} & \\text{{otherwise}} \\end{{cases}}"

        return "?"

    def emit_var(self, name: str) -> str:
        greek = {
            "theta": "\\theta",
            "pi": "\\pi",
            "epsilon": "\\epsilon",
            "delta": "\\delta",
            "PI": "\\pi",
            "EPSILON": "\\epsilon",
        }
        if name in greek:
            return greek[name]
        if len(name) == 1:
            return name
        return f"\\text{{{name}}}"

    def emit_binop(self, op: str, left: str, right: str) -> str:
        if op == "dot":
            return f"{left} \\cdot {right}"
        if op == "cross":
            return f"{left} \\times {right}"
        if op == "/":
            return f"\\frac{{{left}}}{{{right}}}"
        if op == "*":
            return f"{left} \\cdot {right}"
        return f"{left} {op} {right}"

    def emit_call(self, func: str, args: list[str]) -> str:
        if func == "sqrt":
            return f"\\sqrt{{{args[0]}}}"
        if func in ("sin", "cos", "tan"):
            return f"\\{func}({args[0]})"
        if func == "abs":
            return f"|{args[0]}|"
        return f"\\text{{{func}}}({', '.join(args)})"


# ==============================================================================
# Main
# ==============================================================================


def generate_c(funcs: list[FuncDef], mode: Mode, suffix: str = "") -> str:
    """Generate C include file."""
    lines = []
    mode_str = "floating-point" if mode == Mode.FLOAT else "fixed-point"

    lines.append(f"/* Auto-generated from math.dsl - DO NOT EDIT")
    lines.append(f" *")
    lines.append(f" * Math operations for B3D ({mode_str} mode).")
    lines.append(f" * I❤LA-style DSL compiled to C.")
    lines.append(f" */")
    lines.append("")

    gen = CodeGen(mode, suffix)
    for func in funcs:
        lines.append(gen.emit_func(func))
        lines.append("")

    return "\n".join(lines)


def generate_latex(funcs: list[FuncDef]) -> str:
    """Generate LaTeX documentation."""
    lines = []

    lines.append("% Auto-generated from math.dsl")
    lines.append("\\documentclass{article}")
    lines.append("\\usepackage{amsmath}")
    lines.append("\\begin{document}")
    lines.append("")
    lines.append("\\section{B3D Math Functions}")
    lines.append("")
    lines.append("\\begin{align*}")

    gen = LaTeXGen()
    for func in funcs:
        lines.append(gen.emit_func(func))

    lines.append("\\end{align*}")
    lines.append("\\end{document}")

    return "\n".join(lines)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate C/LaTeX from I❤LA-style math DSL"
    )
    parser.add_argument("--dsl", default="src/math.dsl", help="DSL source file")
    parser.add_argument("-o", "--output", help="C output file")
    parser.add_argument("--latex", help="LaTeX output file")
    parser.add_argument("--suffix", default="", help="Function name suffix")
    parser.add_argument(
        "--mode",
        choices=["float", "fixed"],
        default="float",
        help="Numeric mode",
    )
    parser.add_argument("--debug", action="store_true", help="Print debug info")
    args = parser.parse_args()

    with open(args.dsl, "r", encoding="utf-8") as f:
        source = f.read()

    # Tokenize
    lexer = Lexer(source)
    tokens = lexer.tokenize()

    if args.debug:
        print("=== Tokens ===")
        for tok in tokens[:50]:
            print(f"  {tok}")

    # Parse
    parser_obj = Parser(tokens)
    funcs = parser_obj.parse()

    if args.debug:
        print(f"\n=== Parsed {len(funcs)} functions ===")
        for func in funcs:
            print(
                f"  {func.name}: {[p.name for p in func.params]} -> {func.return_type}"
            )

    # Generate C
    if args.output:
        mode = Mode.FLOAT if args.mode == "float" else Mode.FIXED
        output = generate_c(funcs, mode, args.suffix)
        with open(args.output, "w") as f:
            f.write(output)
        print(f"Generated: {args.output}")

    # Generate LaTeX
    if args.latex:
        output = generate_latex(funcs)
        with open(args.latex, "w") as f:
            f.write(output)
        print(f"Generated: {args.latex}")

    # Print to stdout if no output specified
    if not args.output and not args.latex:
        mode = Mode.FLOAT if args.mode == "float" else Mode.FIXED
        print(generate_c(funcs, mode, args.suffix))


if __name__ == "__main__":
    main()
