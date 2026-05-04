"""Parse a small subset of Verilog (assign chains using &, |, ^, ~, parens),
evaluate the resulting Boolean function for all input combinations.
"""
import re
from itertools import product
from typing import Tuple, List, Dict, Any


_TOKEN_RE = re.compile(r'\(|\)|~|&|\||\^|[A-Za-z_][A-Za-z0-9_]*|\s+')


def _tokenize(expr: str) -> List[str]:
    out = []
    for m in _TOKEN_RE.finditer(expr):
        t = m.group(0)
        if t.strip():
            out.append(t)
    return out


def _parse_expr(tokens: List[str]):
    """Recursive-descent parser. Precedence: ~ > & > ^ > |."""
    pos = [0]

    def peek():
        return tokens[pos[0]] if pos[0] < len(tokens) else None

    def consume(expected=None):
        t = peek()
        if expected is not None and t != expected:
            raise ValueError(f"expected {expected!r}, got {t!r}")
        pos[0] += 1
        return t

    def parse_or():
        node = parse_xor()
        while peek() == '|':
            consume('|')
            node = ('or', node, parse_xor())
        return node

    def parse_xor():
        node = parse_and()
        while peek() == '^':
            consume('^')
            node = ('xor', node, parse_and())
        return node

    def parse_and():
        node = parse_unary()
        while peek() == '&':
            consume('&')
            node = ('and', node, parse_unary())
        return node

    def parse_unary():
        if peek() == '~':
            consume('~')
            return ('not', parse_unary())
        return parse_atom()

    def parse_atom():
        t = peek()
        if t == '(':
            consume('(')
            node = parse_or()
            consume(')')
            return node
        if t is None or t in (')', '&', '|', '^', '~'):
            raise ValueError(f"unexpected token: {t!r}")
        consume()
        return ('var', t)

    result = parse_or()
    if pos[0] != len(tokens):
        raise ValueError(f"trailing tokens: {tokens[pos[0]:]}")
    return result


def parse_verilog(path: str) -> dict:
    """Returns {'inputs': [str], 'outputs': [str], 'wires': [str],
                'assigns': [(lhs, expr_ast)]}."""
    with open(path) as f:
        src = f.read()
    src = re.sub(r'//[^\n]*', '', src)
    src = re.sub(r'/\*.*?\*/', '', src, flags=re.DOTALL)

    inputs, outputs, wires, assigns = [], [], [], []

    for stmt in src.split(';'):
        stmt = stmt.strip()
        if not stmt or stmt.startswith('module') or stmt.startswith('endmodule'):
            continue
        if stmt.startswith('input'):
            names = [n.strip() for n in stmt[len('input'):].split(',')]
            inputs.extend(n for n in names if n)
        elif stmt.startswith('output'):
            names = [n.strip() for n in stmt[len('output'):].split(',')]
            outputs.extend(n for n in names if n)
        elif stmt.startswith('wire'):
            names = [n.strip() for n in stmt[len('wire'):].split(',')]
            wires.extend(n for n in names if n)
        elif stmt.startswith('assign'):
            body = stmt[len('assign'):].strip()
            lhs, rhs = body.split('=', 1)
            lhs = lhs.strip()
            tokens = _tokenize(rhs)
            ast = _parse_expr(tokens)
            assigns.append((lhs, ast))
        else:
            raise ValueError(f"unsupported Verilog statement: {stmt[:80]!r}")

    return {'inputs': inputs, 'outputs': outputs,
            'wires': wires, 'assigns': assigns}


def io_signals(path: str) -> Tuple[List[str], List[str]]:
    p = parse_verilog(path)
    return p['inputs'], p['outputs']


def _eval(ast, env: Dict[str, int]) -> int:
    op = ast[0]
    if op == 'var':
        return env[ast[1]]
    if op == 'not':
        return 1 - _eval(ast[1], env)
    a = _eval(ast[1], env)
    b = _eval(ast[2], env)
    if op == 'and':
        return a & b
    if op == 'or':
        return a | b
    if op == 'xor':
        return a ^ b
    raise ValueError(f"unknown op: {op}")


def compute_truth_table(path: str) -> Dict[Tuple[int, ...], Tuple[int, ...]]:
    p = parse_verilog(path)
    inputs = p['inputs']
    outputs = p['outputs']
    table = {}
    for bits in product((0, 1), repeat=len(inputs)):
        env = dict(zip(inputs, bits))
        for lhs, ast in p['assigns']:
            env[lhs] = _eval(ast, env)
        table[bits] = tuple(env[o] for o in outputs)
    return table
