// Fast beam-search solver for CCC level 6.
// Usage: solver2 input output beam threads [wBuilt wMiss wMax wSkirt wJunk]
#include <bits/stdc++.h>
using namespace std;

static const int DIRTO[6][2] = {{2,0},{1,1},{-1,1},{-2,0},{-1,-1},{1,-1}};      // dc,dr: E SE SW W NW NE
static const int PERIM[6][2] = {{1,-1},{2,0},{1,1},{-1,1},{-2,0},{-1,-1}};      // NE E SE SW W NW
static const char OPS[4] = {'-','T','F','B'};

struct Pattern { int br, rows, cols, cnt; vector<string> g; };
struct Weights { double built = 30, miss = 6, mx = 12, skirt = 60, junk = 0.5, skh = 1; };

template<int W> struct BS {
    uint64_t w[W];
    void clear() { memset(w, 0, sizeof w); }
    void set(int i) { w[i >> 6] |= 1ull << (i & 63); }
    bool get(int i) const { return w[i >> 6] >> (i & 63) & 1; }
};

template<int W> struct Solver {
    typedef BS<W> B;
    int R, P, n, cnt, beam;
    Weights wt;
    vector<B> ray[6];
    B pat, skirt, notpat;
    vector<int> patCells, qCells;      // qCells = pattern cells followed by skirting cells
    int K, Q;
    vector<vector<uint8_t>> cost[6];   // cost[t][pos][q]: steps for active track t at pos to cover cell qCells[q]
    static const int MAXQ = 400;
    uint8_t lv1[7][MAXQ], lv2[7][MAXQ], lv3[7][MAXQ];

    Solver(const Pattern& p, int beam, Weights wt) : cnt(p.cnt), beam(beam), wt(wt) {
        R = p.br + 1;
        int cols = 4*R-3, rows = 2*R-1;
        auto inGrid = [&](int c, int r) {
            if (r < 0 || r >= 2*R-1) return false;
            if (((r+c+R+1)%2+2)%2) return false;
            if (r < R) return c+r >= R-1 && c-r <= cols-R+1;
            return c+(2*R-1-r) >= R-1 && c+r-2*R+1 < cols-R+1;
        };
        auto isInterior = [&](int r, int c) {
            if (!inGrid(c, r) || r == 0 || r == rows-1) return false;
            if (r < R) return !(c == R-1-r || c == cols-1-R+1+r);
            return !(c == R-1+r-2*R+2 || c == cols-1+R-1-r);
        };
        map<pair<int,int>,int> id;
        for (int r = 0; r < rows; r++) for (int c = 0; c < cols; c++) if (isInterior(r,c)) { int k = id.size(); id[{r,c}] = k; }
        n = id.size(); assert(n <= 64*W);
        vector<pair<int,int>> perim; int cc = 0, cr = R-1;
        for (int d = 0; d < 6; d++) for (int i = 0; i < R-1; i++) { perim.push_back({cc,cr}); cc += PERIM[d][0]; cr += PERIM[d][1]; }
        P = perim.size();
        for (int t = 0; t < 6; t++) {
            ray[t].resize(P);
            for (int i = 0; i < P; i++) {
                ray[t][i].clear();
                int c = perim[i].first + DIRTO[t][0], r = perim[i].second + DIRTO[t][1];
                while (inGrid(c, r)) { auto it = id.find({r,c}); if (it != id.end()) ray[t][i].set(it->second); c += DIRTO[t][0]; r += DIRTO[t][1]; }
            }
        }
        pat.clear(); skirt.clear();
        for (int r = 0; r < p.rows; r++) for (int c = 0; c < p.cols; c++) {
            auto it = id.find({r,c});
            if (p.g[r][c] == 'M') { pat.set(it->second); patCells.push_back(it->second); }
            else if (p.g[r][c] == 'O' && it != id.end()) {
                for (int d = 0; d < 6; d++) { int rr = r+DIRTO[d][1], c2 = c+DIRTO[d][0];
                    if (rr >= 0 && rr < p.rows && c2 >= 0 && c2 < p.cols && p.g[rr][c2] == 'M') { skirt.set(it->second); break; } }
            }
        }
        for (int i = 0; i < W; i++) notpat.w[i] = ~pat.w[i];
        qCells = patCells; K = patCells.size();
        for (int i = 0; i < n; i++) if (skirt.get(i)) qCells.push_back(i);
        Q = qCells.size(); assert(Q <= MAXQ);
        for (int t = 0; t < 6; t++) {
            cost[t].assign(P, vector<uint8_t>(Q, 200));
            for (int k = 0; k < Q; k++)
                for (int i = 0; i < P; i++) if (ray[t][i].get(qCells[k]))
                    for (int q = 0; q < P; q++) { int d = abs(q-i); d = min(d, P-d); cost[t][q][k] = min<int>(cost[t][q][k], d); }
        }
    }

    struct St { uint8_t pos[6]; uint8_t act; int acc; B m; int parent; uint16_t code; };
    struct Cand { float score; int parent; uint16_t code; };

    // Apply action code to state s -> n (without parent bookkeeping). Returns accepted flag.
    void apply(const St& s, int code, St& n) const {
        memcpy(n.pos, s.pos, 6); n.act = s.act;
        for (int t = 0; t < 6; t++) {
            int o = (code >> (2*t)) & 3;
            if (o == 1) n.act ^= 1 << t;
            else if (o == 2) n.pos[t] = (n.pos[t]+1) % P;
            else if (o == 3) n.pos[t] = (n.pos[t]+P-1) % P;
        }
        B once, twice, thrice; once.clear(); twice.clear(); thrice.clear();
        for (int t = 0; t < 6; t++) if (n.act >> t & 1) {
            const B& L = ray[t][n.pos[t]];
            for (int i = 0; i < W; i++) { thrice.w[i] |= twice.w[i] & L.w[i]; twice.w[i] |= once.w[i] & L.w[i]; once.w[i] |= L.w[i]; }
        }
        bool full = true, clean = true;
        for (int i = 0; i < W; i++) {
            n.m.w[i] = (s.m.w[i] | twice.w[i]) & ~thrice.w[i];
            if (pat.w[i] & ~n.m.w[i]) full = false;
            if (n.m.w[i] & skirt.w[i]) clean = false;
        }
        n.acc = s.acc;
        if (full && clean) { n.acc++; for (int i = 0; i < W; i++) n.m.w[i] &= ~pat.w[i]; }
    }

    void rec(const St& s, int pi, int t, int code, const B& on, const B& tw, const B& th, vector<Cand>& cands) {
        if (t == 6) {
            B m; bool full = true, clean = true;
            int built = 0, sk = 0, junk = 0;
            for (int i = 0; i < W; i++) {
                m.w[i] = (s.m.w[i] | tw.w[i]) & ~th.w[i];
                if (pat.w[i] & ~m.w[i]) full = false;
                uint64_t skm = m.w[i] & skirt.w[i];
                if (skm) clean = false;
                built += __builtin_popcountll(m.w[i] & pat.w[i]);
                sk += __builtin_popcountll(skm);
                junk += __builtin_popcountll(m.w[i] & ~pat.w[i] & ~skirt.w[i]);
            }
            int acc = s.acc;
            const uint8_t *a2 = lv2[6], *a3 = lv3[6];
            double hs = 0, hm = 0, sc;
            if (full && clean) {
                acc++;
                if (acc >= cnt) { cands.push_back({1e30f, pi, (uint16_t)code}); return; }
                for (int k = 0; k < K; k++) { hs += a2[k]; hm = max(hm, (double)a2[k]); }
                sc = acc * 1e5 - junk * wt.junk;
            } else {
                for (int k = 0; k < K; k++) if (!(m.w[qCells[k] >> 6] >> (qCells[k] & 63) & 1)) { hs += a2[k]; hm = max(hm, (double)a2[k]); }
                for (int k = K; k < Q; k++) if (m.w[qCells[k] >> 6] >> (qCells[k] & 63) & 1) { hs += a3[k] * wt.skh; hm = max(hm, (double)a3[k]); }
                sc = acc * 1e5 + built * wt.built - sk * wt.skirt - junk * wt.junk;
            }
            sc -= hs * wt.miss + hm * wt.mx;
            cands.push_back({(float)sc, pi, (uint16_t)code});
            return;
        }
        for (int o = 0; o < 4; o++) {
            int pos = s.pos[t]; bool a = s.act >> t & 1;
            if (o == 1) a = !a; else if (o == 2) pos = (pos+1) % P; else if (o == 3) pos = (pos+P-1) % P;
            const uint8_t* ct = cost[t][pos].data();
            const uint8_t *p1 = lv1[t], *p2 = lv2[t], *p3 = lv3[t];
            uint8_t *b1 = lv1[t+1], *b2 = lv2[t+1], *b3 = lv3[t+1];
            int add = a ? 0 : 1;
            for (int k = 0; k < Q; k++) {
                int v = ct[k] + add, x1 = p1[k], x2 = p2[k], x3 = p3[k];
                if (v < x1) { b3[k] = x2; b2[k] = x1; b1[k] = v; }
                else if (v < x2) { b3[k] = x2; b2[k] = v; b1[k] = x1; }
                else { b3[k] = min(v, x3); b2[k] = x2; b1[k] = x1; }
            }
            if (a) {
                const B& L = ray[t][pos];
                B on2, tw2, th2;
                for (int i = 0; i < W; i++) { th2.w[i] = th.w[i] | (tw.w[i] & L.w[i]); tw2.w[i] = tw.w[i] | (on.w[i] & L.w[i]); on2.w[i] = on.w[i] | L.w[i]; }
                rec(s, pi, t+1, code | o << (2*t), on2, tw2, th2, cands);
            } else rec(s, pi, t+1, code | o << (2*t), on, tw, th, cands);
        }
    }

    int maxDepth = 150;
    vector<array<char,6>> solve() {
        vector<vector<St>> layers;
        St s0; s0.act = 0; s0.acc = 0; s0.parent = -1; s0.code = 0; s0.m.clear();
        for (int t = 0; t < 6; t++) s0.pos[t] = t * (R-1);
        layers.push_back({s0});
        vector<Cand> cands;
        for (int depth = 0; depth < maxDepth; depth++) {
            const vector<St>& cur = layers.back();
            cands.clear();
            for (int pi = 0; pi < (int)cur.size(); pi++) {
                const St& s = cur[pi];
                B o0, t0, h0; o0.clear(); t0.clear(); h0.clear();
                for (int k = 0; k < Q; k++) lv1[0][k] = lv2[0][k] = lv3[0][k] = 250;
                rec(s, pi, 0, 0, o0, t0, h0, cands);
            }
            // finished?
            for (auto& c : cands) if (c.score > 1e29f) {
                St n; apply(layers.back()[c.parent], c.code, n); n.parent = c.parent; n.code = c.code;
                vector<array<char,6>> prog;
                int code = c.code, p = c.parent;
                for (int d = (int)layers.size(); d >= 1; d--) {
                    array<char,6> l; for (int t = 0; t < 6; t++) l[t] = OPS[(code >> (2*t)) & 3];
                    prog.push_back(l);
                    if (d >= 2) { code = layers[d-1][p].code; p = layers[d-1][p].parent; }
                }
                reverse(prog.begin(), prog.end());
                return prog;
            }
            size_t keep = min(cands.size(), (size_t)beam * 3);
            nth_element(cands.begin(), cands.begin() + keep - 1, cands.end(), [](const Cand& a, const Cand& b){ return a.score > b.score; });
            cands.resize(keep);
            sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b){ return a.score > b.score; });
            vector<St> next; next.reserve(beam);
            unordered_set<uint64_t> seen;
            for (auto& c : cands) {
                St n; apply(cur[c.parent], c.code, n); n.parent = c.parent; n.code = c.code;
                uint64_t h = 1469598103934665603ull;
                for (int i = 0; i < W; i++) h = (h ^ n.m.w[i]) * 1099511628211ull;
                for (int t = 0; t < 6; t++) h = (h ^ n.pos[t]) * 1099511628211ull;
                h = (h ^ n.act) * 1099511628211ull; h = (h ^ n.acc) * 1099511628211ull;
                if (!seen.insert(h).second) continue;
                next.push_back(n);
                if ((int)next.size() >= beam) break;
            }
            layers.push_back(move(next));
        }
        return {};
    }
};

int gMaxDepth = 150;
vector<array<char,6>> solveOne(const Pattern& p, int beam, Weights wt) {
    int R = p.br + 1, n = 3*(R-1)*(R-2)+1;
    if (n <= 128) { Solver<2> s(p, beam, wt); s.maxDepth = gMaxDepth; return s.solve(); }
    if (n <= 256) { Solver<4> s(p, beam, wt); s.maxDepth = gMaxDepth; return s.solve(); }
    Solver<8> s(p, beam, wt); s.maxDepth = gMaxDepth; return s.solve();
}

int main(int argc, char** argv) {
    ifstream in(argv[1]);
    int beam = atoi(argv[3]), threads = atoi(argv[4]);
    Weights wt;
    if (argc > 9) { wt.built = atof(argv[5]); wt.miss = atof(argv[6]); wt.mx = atof(argv[7]); wt.skirt = atof(argv[8]); wt.junk = atof(argv[9]); }
    if (argc > 10) wt.skh = atof(argv[10]);
    int N; in >> N;
    vector<Pattern> pats(N);
    for (auto& p : pats) { in >> p.br >> p.rows >> p.cols >> p.cnt; p.g.resize(p.rows); for (auto& l : p.g) in >> l; }
    vector<vector<array<char,6>>> res(N);
    vector<int> order;
    const char* sel = getenv("SEL");      // comma-separated pattern indices
    if (sel) { stringstream ss(sel); string x; while (getline(ss, x, ',')) order.push_back(stoi(x)); }
    else for (int i = 0; i < N; i++) order.push_back(i);
    if (getenv("MAXD")) gMaxDepth = atoi(getenv("MAXD"));
    const char* dir = getenv("SOLDIR");
    atomic<int> next{0};
    vector<thread> th;
    for (int k = 0; k < threads; k++) th.emplace_back([&]{
        int j;
        while ((j = next++) < (int)order.size()) {
            int i = order[j];
            res[i] = solveOne(pats[i], beam, wt);
            fprintf(stderr, "pattern %d cnt=%d lines=%zu\n", i, pats[i].cnt, res[i].size());
            if (dir && !res[i].empty()) {
                string f = string(dir) + "/" + to_string(i) + "_" + to_string(res[i].size()) + "_" + to_string(getpid()) + ".txt";
                ofstream o(f); for (auto& l : res[i]) for (int t = 0; t < 6; t++) o << l[t] << (t < 5 ? ' ' : '\n');
            }
        }
    });
    for (auto& t : th) t.join();
    ofstream out(argv[2]);
    long total = 0;
    for (int i = 0; i < N; i++) {
        if (res[i].empty()) res[i].push_back({'-','-','-','-','-','-'});
        for (auto& l : res[i]) for (int t = 0; t < 6; t++) out << l[t] << (t < 5 ? ' ' : '\n');
        if (i + 1 < N) out << '\n';
        total += res[i].size();
    }
    fprintf(stderr, "TOTAL %ld\n", total);
}
