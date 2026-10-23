"""Faithful port of the level-6 sandbox simulator. Usage: verify.py input output"""
import sys

DIR_TO = [(2, 0), (1, 1), (-1, 1), (-2, 0), (-1, -1), (1, -1)]  # (dc, dr)
PERIM_DIRS = [(1, -1), (2, 0), (1, 1), (-1, 1), (-2, 0), (-1, -1)]  # NE E SE SW W NW
FIRE_OV = DIR_TO  # E SE SW W NW NE


def parse_input(text):
    lines = text.split('\n')
    n = int(lines[0])
    pats, i = [], 1
    for _ in range(n):
        while not lines[i].strip():
            i += 1
        br, rows, cols, cnt = map(int, lines[i].split())
        grid = [list(l) for l in lines[i + 1:i + 1 + rows]]
        pats.append((br, rows, cols, cnt, grid))
        i += 1 + rows
    return pats


def parse_output(text, n):
    progs = [p for p in text.strip('\n').split('\n\n')]
    progs = [[l.split() for l in p.split('\n') if l.strip()] for p in progs]
    assert len(progs) == n, (len(progs), n)
    return progs


def simulate(pat, prog):
    br, prow, pcol, cnt, pgrid = pat
    R = br + 1
    cols, rows = 4 * R - 3, 2 * R - 1
    perim, cur = [], (0, R - 1)
    for dc, dr in PERIM_DIRS:
        for _ in range(R - 1):
            perim.append(cur)
            cur = (cur[0] + dc, cur[1] + dr)

    def in_grid(c, r):
        if r < 0 or r >= 2 * R - 1: return False
        if (r + c + R + 1) % 2: return False
        if r < R: return c + r >= R - 1 and c - r <= cols - R + 1
        return c + (2 * R - 1 - r) >= R - 1 and c + r - 2 * R + 1 < cols - R + 1

    def base(r, c):
        if (r + c + R + 1) % 2: return '-'
        rr = r == 0 or r == rows - 1
        if r < R:
            cr = c == R - 1 - r or c == cols - 1 - R + 1 + r
            if c + r >= R - 1 and c - r <= cols - R + 1: return '#' if rr or cr else 'O'
            return ':'
        cr = c == R - 1 + r - 2 * R + 2 or c == cols - 1 + R - 1 - r
        if c + 2 * R - 1 - r >= R - 1 and c + r - 2 * R + 1 < cols - R + 1: return '#' if rr or cr else 'O'
        return ':'

    grid = [[base(r, c) for c in range(cols)] for r in range(rows)]
    ray = []
    for (c0, r0) in perim:
        per = []
        for dc, dr in FIRE_OV:
            cells, c, r = [], c0 + dc, r0 + dr
            while in_grid(c, r):
                cells.append((r, c)); c += dc; r += dr
            per.append(cells)
        ray.append(per)
    P = len(perim)
    idx = [i * (R - 1) for i in range(6)]
    act = [False] * 6
    acc = 0

    def skirting(r, c):
        return any(0 <= r + dr < prow and 0 <= c + dc < pcol and pgrid[r + dr][c + dc] == 'M' for dc, dr in DIR_TO)

    for line in prog:
        for r in range(rows):
            for c in range(cols):
                if grid[r][c] in 'xd': grid[r][c] = 'O'
                elif grid[r][c] == 'm': grid[r][c] = 'M'
        line = (line + ['-'] * 6)[:6]
        for t in range(6):
            if line[t] == 'F': idx[t] = (idx[t] + 1) % P
            elif line[t] == 'B': idx[t] = (idx[t] + P - 1) % P
            elif line[t] == 'T': act[t] = not act[t]
            else: assert line[t] == '-'
        hit = {}
        for t in range(6):
            if act[t]:
                for k in ray[idx[t]][t]: hit[k] = hit.get(k, 0) + 1
        for (r, c), n in hit.items():
            cell = grid[r][c]
            if cell not in 'OM': continue
            if n == 2: grid[r][c] = 'm' if cell == 'O' else 'M'
            elif n >= 3: grid[r][c] = 'd'
        ok = True
        for r in range(prow):
            for c in range(pcol):
                p, g = pgrid[r][c], grid[r][c]
                if p == 'M' and g not in 'Mm': ok = False
                elif p == 'O' and skirting(r, c) and g in 'Mm': ok = False
        if ok:
            for r in range(prow):
                for c in range(pcol):
                    if pgrid[r][c] == 'M' and grid[r][c] in 'Mm': grid[r][c] = 'x'
            acc += 1
    return acc


def main():
    pats = parse_input(open(sys.argv[1]).read())
    progs = parse_output(open(sys.argv[2]).read(), len(pats))
    total, bad = 0, 0
    for i, (pat, prog) in enumerate(zip(pats, progs)):
        acc = simulate(pat, prog)
        total += len(prog)
        if acc < pat[3]:
            bad += 1
            print(f'pattern {i}: accepted {acc}/{pat[3]}')
    print(f'lines={total} failed={bad}')


if __name__ == '__main__':
    main()
