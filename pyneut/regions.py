"""Region expressions: the cell-geometry syntax of the input format.

    "-shell"                     inside surface `shell`
    "+floor"                     outside surface `floor`
    "-outer & -top & +bottom"    intersection
    "-a | -b"                    union
    "~(-hole)"                   complement
    "-ball & ~(-hole | -slot)"   parentheses group

Precedence: ``~`` binds tightest, then ``&``, then ``|``. Intersection must be
written with ``&``: a bare space between two terms is an error, so a missing
operator is never silently read as an intersection.

An expression is parsed into a small tree (`parse_region`) and then compiled
onto the existing CSG objects (`compile_region`): a ``-s`` half-space is the
surface itself, ``+s`` and ``~x`` are complement Regions, ``&`` and ``|`` are
intersection and union Regions. The transport kernels are unchanged.
"""
import re
from dataclasses import dataclass

from .medium import Region

_NAME = r"[A-Za-z_][A-Za-z0-9_]*"
_TOKEN_RE = re.compile(rf"\s*(?:(?P<half>[-+]{_NAME})|(?P<op>[&|~()]))")


@dataclass(frozen=True)
class HalfSpace:
    surface: str
    inside: bool          # True for "-name", False for "+name"


@dataclass(frozen=True)
class And:
    terms: tuple


@dataclass(frozen=True)
class Or:
    terms: tuple


@dataclass(frozen=True)
class Not:
    term: object


def _tokenize(text):
    tokens, pos = [], 0
    text = text.rstrip()
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if not m:
            bad = text[pos:].strip().split()[0] if text[pos:].strip() else text[pos:]
            hint = ""
            if re.fullmatch(_NAME, bad):
                hint = f" (a surface needs a sign: '-{bad}' for inside, '+{bad}' for outside)"
            raise ValueError(f"unexpected {bad!r} at position {pos + 1}{hint}")
        tokens.append(m.group("half") or m.group("op"))
        pos = m.end()
    return tokens


class _Parser:
    def __init__(self, tokens, text):
        self.tokens, self.i, self.text = tokens, 0, text

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def parse(self):
        if not self.tokens:
            raise ValueError("empty region expression")
        node = self.union()
        tok = self.peek()
        if tok is not None:
            if tok == ")":
                raise ValueError("unbalanced ')'")
            prev = self.tokens[self.i - 1]
            raise ValueError(f"missing operator between {prev!r} and {tok!r}; "
                             "write '&' for intersection or '|' for union")
        return node

    def union(self):
        terms = [self.intersection()]
        while self.peek() == "|":
            self.take()
            terms.append(self.intersection())
        return terms[0] if len(terms) == 1 else Or(tuple(terms))

    def intersection(self):
        terms = [self.factor()]
        while self.peek() == "&":
            self.take()
            terms.append(self.factor())
        return terms[0] if len(terms) == 1 else And(tuple(terms))

    def factor(self):
        tok = self.take()
        if tok is None:
            raise ValueError("expression ends where a term was expected")
        if tok == "~":
            return Not(self.factor())
        if tok == "(":
            node = self.union()
            if self.take() != ")":
                raise ValueError("missing ')'")
            return node
        if tok in ("&", "|", ")"):
            raise ValueError(f"unexpected {tok!r} where a term was expected")
        return HalfSpace(tok[1:], tok[0] == "-")


def parse_region(text):
    """Parse a region expression into a tree of HalfSpace/And/Or/Not."""
    if not isinstance(text, str):
        raise ValueError(f"a region is a string expression, got {text!r}")
    return _Parser(_tokenize(text), text).parse()


def surface_names(node):
    """Every surface name an expression tree refers to."""
    if isinstance(node, HalfSpace):
        return {node.surface}
    if isinstance(node, Not):
        return surface_names(node.term)
    return set().union(*(surface_names(t) for t in node.terms))


def compile_region(node, surfaces):
    """Compile an expression tree onto CSG objects.

    `surfaces` maps names to primitive surfaces (Plane, Sphere, Cylinder) or
    Box regions. Returns a surface or a Region usable inside a Region's
    surface list.
    """
    if isinstance(node, HalfSpace):
        surface = surfaces[node.surface]
        return surface if node.inside else Region(surfaces=[surface],
                                                  operation="complement")
    if isinstance(node, Not):
        return Region(surfaces=[compile_region(node.term, surfaces)],
                      operation="complement")
    op = "intersection" if isinstance(node, And) else "union"
    return Region(surfaces=[compile_region(t, surfaces) for t in node.terms],
                  operation=op)
