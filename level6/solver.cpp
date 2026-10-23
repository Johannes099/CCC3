// Beam-search solver for CCC level 6. Usage: solver input output [beam] [threads]
#include <bits/stdc++.h>
using namespace std;
typedef bitset<512> BS;

static const int DIRTO[6][2] = {{2,0},{1,1},{-1,1},{-2,0},{-1,-1},{1,-1}};      // dc,dr: E SE SW W NW NE
static const int PERIM[6][2] = {{1,-1},{2,0},{1,1},{-1,1},{-2,0},{-1,-1}};      // NE E SE SW W NW

struct Pattern { int br, rows, cols, cnt; vector<string> g; };

struct Board {
    int R, P, n;                       // radius, perimeter length, interior cell count
    vector<BS> ray[6];                 // ray[t][idx]
    BS pat, skirt;                     // pattern matter cells, skirting cells
    vector<int> patCells;
    vector<vector<int>> need[6];       // need[t][k] = perimeter indices whose ray from t covers patCells[k]
};

Board build(const Pattern& p) {
    Board b; int R = p.br + 1; b.R = R;
    int cols = 4*R-3, rows = 2*R-1;
    auto inGrid = [&](int c, int r) {
        if (r < 0 || r >= 2*R-1) return false;
        if (((r+c+R+1)%2+2)%2) return false;
        if (r < R) return c+r >= R-1 && c-r <= cols-R+1;
        return c+(2*R-1-r) >= R-1 && c+r-2*R+1 < cols-R+1;
    };
    auto isInterior = [&](int r, int c) {
        if (!inGrid(c, r)) return false;
        if (r == 0 || r == rows-1) return false;
        if (r < R) return !(c == R-1-r || c == cols-1-R+1+r);
        return !(c == R-1+r-2*R+2 || c == cols-1+R-1-r);
    };
    map<pair<int,int>,int> id;
    for (int r = 0; r < rows; r++) for (int c = 0; c < cols; c++) if (isInterior(r,c)) id[{r,c}] = (int)id.size();
    b.n = id.size();
    vector<pair<int,int>> perim; int cc = 0, cr = R-1;
    for (int d = 0; d < 6; d++) for (int i = 0; i < R-1; i++) { perim.push_back({cc,cr}); cc += PERIM[d][0]; cr += PERIM[d][1]; }
    b.P = perim.size();
    for (int t = 0; t < 6; t++) {
        b.ray[t].resize(b.P);
        for (int i = 0; i < b.P; i++) {
            int c = perim[i].first + DIRTO[t][0], r = perim[i].second + DIRTO[t][1];
            while (inGrid(c, r)) { auto it = id.find({r,c}); if (it != id.end()) b.ray[t][i][it->second] = 1; c += DIRTO[t][0]; r += DIRTO[t][1]; }
        }
    }
    for (int r = 0; r < p.rows; r++) for (int c = 0; c < p.cols; c++) {
        auto it = id.find({r,c});
        if (p.g[r][c] == 'M') { assert(it != id.end()); b.pat[it->second] = 1; b.patCells.push_back(it->second); }
        else if (p.g[r][c] == 'O' && it != id.end()) {
            for (int d = 0; d < 6; d++) { int rr = r+DIRTO[d][1], c2 = c+DIRTO[d][0];
                if (rr >= 0 && rr < p.rows && c2 >= 0 && c2 < p.cols && p.g[rr][c2] == 'M') { b.skirt[it->second] = 1; break; } }
        }
    }
    for (int t = 0; t < 6; t++) {
        b.need[t].resize(b.patCells.size());
        for (size_t k = 0; k < b.patCells.size(); k++)
            for (int i = 0; i < b.P; i++) if (b.ray[t][i][b.patCells[k]]) b.need[t][k].push_back(i);
    }
    return b;
}

struct State { array<uint8_t,6> pos; uint8_t act; BS m; int acc; double score; int parent; array<char,6> op; };

struct Solver {
    const Board& b; int cnt; int beam;
    // weights
    double wAcc = 10000, wBuilt = 30, wMiss = 6, wSkirt = 60, wJunk = 0.5;
    Solver(const Board& b, int cnt, int beam) : b(b), cnt(cnt), beam(beam) {}

    int cdist(int a, int c) const { int d = abs(a-c) % b.P; return min(d, b.P-d); }

    double evaluate(const State& s) const {
        double sc = s.acc * wAcc;
        BS miss = b.pat & ~s.m;
        int built = (b.pat & s.m).count();
        sc += built * wBuilt;
        sc -= (s.m & b.skirt).count() * wSkirt;
        sc -= (s.m & ~b.pat & ~b.skirt).count() * wJunk;
        double hsum = 0, hmax = 0;
        for (size_t k = 0; k < b.patCells.size(); k++) {
            if (!miss[b.patCells[k]]) continue;
            int c1 = 1000, c2 = 1000;
            for (int t = 0; t < 6; t++) {
                int best = 1000;
                for (int i : b.need[t][k]) best = min(best, cdist(s.pos[t], i));
                if (best == 1000) continue;
                best += ((s.act >> t) & 1) ? 0 : 1;
                if (best < c1) { c2 = c1; c1 = best; } else if (best < c2) c2 = best;
            }
            hsum += c2; hmax = max(hmax, (double)c2);
        }
        sc -= hsum * wMiss + hmax * wMiss * 2;
        return sc;
    }

    vector<array<char,6>> solve() {
        vector<vector<State>> layers;
        State s0; s0.act = 0; s0.acc = 0; s0.parent = -1;
        for (int t = 0; t < 6; t++) s0.pos[t] = t * (b.R-1);
        s0.score = evaluate(s0);
        layers.push_back({s0});
        const char OPS[4] = {'-','T','F','B'};
        for (int depth = 0; depth < 400; depth++) {
            const vector<State>& cur = layers.back();
            vector<State> next; next.reserve(cur.size() * 64);
            unordered_set<size_t> seen;
            for (int pi = 0; pi < (int)cur.size(); pi++) {
                const State& s = cur[pi];
                for (int code = 0; code < 4096; code++) {
                    State n; n.pos = s.pos; n.act = s.act; n.parent = pi;
                    int cc = code; bool useless = false;
                    for (int t = 0; t < 6; t++) {
                        int o = cc & 3; cc >>= 2; n.op[t] = OPS[o];
                        if (o == 1) n.act ^= 1 << t;
                        else if (o == 2) n.pos[t] = (n.pos[t]+1) % b.P;
                        else if (o == 3) n.pos[t] = (n.pos[t]+b.P-1) % b.P;
                    }
                    if (useless) continue;
                    BS once, twice, thrice;
                    for (int t = 0; t < 6; t++) if ((n.act >> t) & 1) {
                        const BS& L = b.ray[t][n.pos[t]];
                        thrice |= twice & L; twice |= once & L; once |= L;
                    }
                    BS ex2 = twice & ~thrice;
                    n.m = (s.m | ex2) & ~thrice;
                    n.acc = s.acc;
                    if ((b.pat & ~n.m).none() && (n.m & b.skirt).none()) { n.acc++; n.m &= ~b.pat; }
                    size_t h = hash<BS>()(n.m) * 1315423911u ^ n.act * 2654435761u ^ (size_t)n.acc * 97531;
                    for (int t = 0; t < 6; t++) h = h * 131 + n.pos[t];
                    if (!seen.insert(h).second) continue;
                    if (n.acc >= cnt) {
                        // reconstruct
                        vector<array<char,6>> prog{n.op};
                        int p = n.parent;
                        for (int d = (int)layers.size()-1; d >= 1; d--) { prog.push_back(layers[d][p].op); p = layers[d][p].parent; }
                        reverse(prog.begin(), prog.end());
                        return prog;
                    }
                    n.score = evaluate(n);
                    next.push_back(n);
                }
            }
            if ((int)next.size() > beam) {
                nth_element(next.begin(), next.begin()+beam, next.end(), [](const State& a, const State& c){ return a.score > c.score; });
                next.resize(beam);
            }
            layers.push_back(move(next));
        }
        return {};
    }
};

int main(int argc, char** argv) {
    ifstream in(argv[1]);
    int beam = argc > 3 ? atoi(argv[3]) : 50;
    int threads = argc > 4 ? atoi(argv[4]) : 4;
    int N; in >> N;
    vector<Pattern> pats(N);
    for (auto& p : pats) { in >> p.br >> p.rows >> p.cols >> p.cnt; p.g.resize(p.rows); for (auto& l : p.g) in >> l; }
    vector<vector<array<char,6>>> res(N);
    atomic<int> next{0};
    vector<thread> th;
    for (int k = 0; k < threads; k++) th.emplace_back([&]{
        int i;
        while ((i = next++) < N) {
            Board b = build(pats[i]);
            Solver s(b, pats[i].cnt, beam);
            res[i] = s.solve();
            fprintf(stderr, "pattern %d R=%d cnt=%d cells=%zu lines=%zu\n", i, b.R, pats[i].cnt, b.patCells.size(), res[i].size());
        }
    });
    for (auto& t : th) t.join();
    ofstream out(argv[2]);
    long total = 0;
    for (int i = 0; i < N; i++) {
        for (auto& l : res[i]) { for (int t = 0; t < 6; t++) out << l[t] << (t < 5 ? ' ' : '\n'); }
        if (i + 1 < N) out << '\n';
        total += res[i].size();
    }
    fprintf(stderr, "TOTAL %ld\n", total);
}
