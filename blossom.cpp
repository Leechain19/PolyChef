//
// Created by AnthonyZhang on 2025/1/6.
//

#include "blossom.h"
#include <numeric>

int blossom::Blossom(std::vector<std::vector<int>>& graph) {
    int n = static_cast<int>(graph.size()), timer = -1;
    std::vector<int> match(n, -1), label(n), parent(n), orig(n), aux(n, -1), q;
    auto lca = [&](int x, int y) {
        for (timer++; ; std::swap(x, y)) {
            if (x == -1) continue;
            if (aux[x] == timer) return x;
            aux[x] = timer;
            x = (match[x] == -1 ? -1 : orig[parent[match[x]]]);
        }
    };
    auto blossom = [&](int v, int w, int a) {
        while (orig[v] != a) {
            parent[v] = w; w = match[v];
            if (label[w] == 1) label[w] = 0, q.push_back(w);
            orig[v] = orig[w] = a; v = parent[w];
        }
    };
    auto augment = [&](int v) {
        while (v != -1) {
            int pv = parent[v], nv = match[pv];
            match[v] = pv; match[pv] = v; v = nv;
        }
    };
    auto bfs = [&](int root) {
        fill(label.begin(), label.end(), -1);
        std::iota(orig.begin(), orig.end(), 0);
        q.clear();
        label[root] = 0; q.push_back(root);
        for (int i = 0; i < (int)q.size(); ++i) {
            int v = q[i];
            for (auto x : graph[v]) {
                if (label[x] == -1) {
                    label[x] = 1; parent[x] = v;
                    if (match[x] == -1) return augment(x), 1;
                    label[match[x]] = 0; q.push_back(match[x]);
                }
                else if (label[x] == 0 && orig[v] != orig[x]) {
                    int a = lca(orig[v], orig[x]);
                    blossom(x, v, a); blossom(v, x, a);
                }
            }
        }
        return 0;
    };

    for (int i = 0; i < n; i ++) {
        if (match[i] == -1) {
            bfs(i);
        }
    }

    int ret = 0;
    for (auto x : match) {
        ret += (x >= 0);
    }
    return ret;
}