# -*- coding: utf-8 -*-
"""
Бинарная последовательная композиция конечных автоматов A∙B и мутационные тесты
для хвостовой компоненты B.

Модель: внешние входы I подаются только на A, выход A (внутренний канал U)
является входом B, наблюдается только выход B (канал O). Автомат задается
парой таблиц (delta, lam): delta[s][i] - следующее состояние, lam[s][i] - выход;
начальное состояние 0. Сигнал сброса возвращает композицию в начальное
состояние; длина теста считается с учетом сброса (см. length).

Мутант компоненты B задается списком неисправностей (b, u, вид, значение):
вид 'T' - ошибка перехода (переход (b, u) ведет в состояние значение),
вид 'O' - ошибка выхода (на переходе (b, u) выдается выход значение).
Порядок мутанта - число неисправностей в списке.

Тесты:
  обход (cpp_tour)     - кратчайший обход всех переходов автомата Spec = A∙B
                         (задача китайского почтальона, сброс = 1 символ);
  покрытие (comp_cover) - жадное покрытие переходов B, используемых в композиции;
  mb()                 - дополнение теста кратчайшими различающими
                         последовательностями для необнаруженных мутантов
                         (MB1 = обход + мутанты первого порядка,
                          MB12 = MB1 + мутанты второго порядка,
                          MBC = покрытие + мутанты первого порядка).
"""
import itertools
from collections import deque, defaultdict

import networkx as nx


def rand_fsm(n, nin, nout, rng):
    """Случайный полностью определенный детерминированный автомат:
    n состояний, nin входов, nout выходов (равномерный независимый выбор)."""
    d = [[rng.randrange(n) for _ in range(nin)] for _ in range(n)]
    l = [[rng.randrange(nout) for _ in range(nin)] for _ in range(n)]
    return d, l


class Net:
    """Композиция A∙B; хранит достижимую часть автомата Spec = A∙B."""

    def __init__(self, A, B, p, pU=None, pO=None):
        self.dA, self.lA = A
        self.dB, self.lB = B
        self.p = p                                   # |I|
        self.pU = pU if pU is not None else p        # |U|
        self.pO = pO if pO is not None else p        # |O|
        self.N, self.K = len(self.dA), len(self.dB)
        idx = {(0, 0): 0}
        st = [(0, 0)]
        q = deque([(0, 0)])
        nxt, lab = [], []
        while q:
            a, b = q.popleft()
            row_n, row_l = [], []
            for i in range(p):
                u = self.lA[a][i]
                t = (self.dA[a][i], self.dB[b][u])
                if t not in idx:
                    idx[t] = len(st); st.append(t); q.append(t)
                row_n.append(idx[t]); row_l.append((b, u))
            nxt.append(row_n); lab.append(row_l)
        self.st, self.idx, self.nxt, self.lab = st, idx, nxt, lab
        self.r = len(st)                                        # достижимых состояний Spec
        self.used = set(l for row in lab for l in row)          # переходы B, используемые в композиции
        self.reachB = set(b for (_, b) in st)                   # состояния B, достижимые в композиции

    # ---------- мутанты ----------
    def mutant(self, faults):
        dB = [row[:] for row in self.dB]
        lB = [row[:] for row in self.lB]
        for (b, u, k, v) in faults:
            if k == 'T': dB[b][u] = v
            else: lB[b][u] = v
        return dB, lB

    def first_order(self, kinds=('T',)):
        """все мутанты первого порядка заданных видов (списки из одной неисправности)"""
        res = []
        for b in range(self.K):
            for u in range(self.pU):
                if 'T' in kinds:
                    for t in range(self.K):
                        if t != self.dB[b][u]: res.append(((b, u, 'T', t),))
                if 'O' in kinds:
                    for o in range(self.pO):
                        if o != self.lB[b][u]: res.append(((b, u, 'O', o),))
        return res

    # ---------- моделирование ----------
    def run(self, seq, M):
        """номер (с 1) первого символа, на котором выходы Spec и A∙M различаются; 0 - не различаются"""
        dA, lA, dB, lB = self.dA, self.lA, self.dB, self.lB
        dM, lM = M
        a = b = c = 0
        for k, i in enumerate(seq):
            u = lA[a][i]
            if lB[b][u] != lM[c][u]:
                return k + 1
            a, b, c = dA[a][i], dB[b][u], dM[c][u]
        return 0

    def detects(self, test, M):
        return any(self.run(s, M) for s in test)

    def end_triple(self, seq, M):
        dA, lA, dB = self.dA, self.lA, self.dB
        dM = M[0]
        a = b = c = 0
        for i in seq:
            u = lA[a][i]
            a, b, c = dA[a][i], dB[b][u], dM[c][u]
        return (a, b, c)

    def shortest_dist(self, test, M, ext=True):
        """Кратчайший способ обнаружить мутант M: поиск в ширину по произведению Spec и A∙M.
        Возвращает (стоимость, j, суффикс): j = -1 - новая последовательность после сброса
        (стоимость = длина + 1), j >= 0 - продолжение последовательности j теста (только при ext=True).
        None, если A∙M эквивалентен Spec (мутант не обнаружим через контекст)."""
        dA, lA, dB, lB = self.dA, self.lA, self.dB, self.lB
        dM, lM = M
        dist, par = {}, {}
        buckets = defaultdict(list)
        if ext:
            for j, s in enumerate(test):
                t = self.end_triple(s, M)
                if t not in dist:
                    dist[t] = 0; par[t] = ('leaf', j); buckets[0].append(t)
        t0 = (0, 0, 0)
        if t0 not in dist:
            dist[t0] = 1; par[t0] = ('root', -1); buckets[1].append(t0)
        d = 0
        maxd = 1
        while d <= maxd:
            for t in buckets.get(d, ()):
                if dist[t] != d: continue
                a, b, c = t
                for i in range(self.p):
                    u = lA[a][i]
                    if lB[b][u] != lM[c][u]:
                        suf = [i]; x = t
                        while par[x][0] == 'e':
                            suf.append(par[x][2]); x = par[x][1]
                        suf.reverse()
                        kind, j = par[x]
                        return (d + 1, j, suf)
                    nt = (dA[a][i], dB[b][u], dM[c][u])
                    if nt not in dist or dist[nt] > d + 1:
                        dist[nt] = d + 1; par[nt] = ('e', t, i)
                        buckets[d + 1].append(nt); maxd = max(maxd, d + 1)
            d += 1
        return None

    def distinguishable(self, M):
        return self.shortest_dist([], M, False) is not None

    # ---------- обход и покрытие ----------
    def cpp_tour(self):
        """Кратчайший обход всех переходов достижимой части Spec со сбросами (задача китайского почтальона)."""
        r, p = self.r, self.p
        R = r
        bal = [0] * (r + 1)
        for s in range(r):
            for i in range(p):
                t = self.nxt[s][i]
                bal[s] += 1; bal[t] -= 1
        bal[R] += 1; bal[0] -= 1
        G = nx.DiGraph()
        for v in range(r + 1):
            G.add_node(v, demand=bal[v])
        for s in range(r):
            for i in range(p):
                t = self.nxt[s][i]
                if t != s and not G.has_edge(s, t):
                    G.add_edge(s, t, weight=1000, inp=i)
            G.add_edge(s, R, weight=0)
        G.add_edge(R, 0, weight=1001)
        flow = nx.min_cost_flow(G)
        MG = nx.MultiDiGraph()
        for s in range(r):
            for i in range(p):
                MG.add_edge(s, self.nxt[s][i], inp=i)
        MG.add_edge(R, 0, inp=None)
        for u, fl in flow.items():
            for v, f in fl.items():
                for _ in range(f):
                    MG.add_edge(u, v, inp=G[u][v].get('inp'))
        seqs, cur = [], None
        for (u, v, k) in nx.eulerian_circuit(MG, source=R, keys=True):
            if u == R:
                cur = []
            elif v == R:
                seqs.append(cur); cur = None
            else:
                cur.append(MG[u][v][k]['inp'])
        return seqs

    def _nearest(self, src, unc):
        seen = {src: None}
        q = deque([src])
        while q:
            s = q.popleft()
            for i in range(self.p):
                if self.lab[s][i] in unc:
                    path = [i]; x = s
                    while seen[x] is not None:
                        x, ii = seen[x]; path.append(ii)
                    path.reverse()
                    return path
            for i in range(self.p):
                t = self.nxt[s][i]
                if t not in seen:
                    seen[t] = (s, i); q.append(t)
        return None

    def comp_cover(self):
        """Жадное покрытие всех переходов B, используемых в композиции."""
        unc = set(self.used)
        seqs, cur, cs = [], None, 0
        while unc:
            best = None
            if cur is not None:
                pth = self._nearest(cs, unc)
                if pth is not None: best = (len(pth), 'cont', pth)
            pth0 = self._nearest(0, unc)
            if best is None or len(pth0) + 1 < best[0]:
                best = (len(pth0) + 1, 'new', pth0)
            _, kind, pth = best
            if kind == 'new':
                cur = []; seqs.append(cur); cs = 0
            for i in pth:
                unc.discard(self.lab[cs][i]); cur.append(i); cs = self.nxt[cs][i]
        return seqs


def length(test):
    """длина теста: сумма длин последовательностей плюс один сброс на каждую"""
    return sum(len(s) + 1 for s in test)


def mb(net, init, mutants, ext):
    """Мутационный тест: к начальному тесту init добавляются кратчайшие различающие
    последовательности для мутантов, которые тест еще не обнаруживает
    (ext=False - новая последовательность после сброса; ext=True - можно продолжить существующую)."""
    test = [list(s) for s in init]
    added = 0
    for M in mutants:
        if net.detects(test, M): continue
        r = net.shortest_dist(test, M, ext)
        if r is None: continue
        cost, j, suf = r
        if j < 0: test.append(suf)
        else: test[j].extend(suf)
        added += 1
    return test, added


# ---------- трансляция через головную компоненту ----------
def translate(A, beta, p):
    """Трансляция последовательности beta (в алфавите U) через автомат A (обратный автомат):
    входная последовательность A, порождающая наибольший префикс beta, который A может выдать
    из начального состояния."""
    dA, lA = A
    cur = {0: []}
    for u in beta:
        nxt = {}
        for a, alpha in cur.items():
            for i in range(p):
                if lA[a][i] == u and dA[a][i] not in nxt:
                    nxt[dA[a][i]] = alpha + [i]
        if not nxt: break
        cur = nxt
    return next(iter(cur.values()))


# ---------- выборки кратных мутантов (ошибки переходов) ----------
def second_order_sample(net, rng, size):
    """Мутанты второго порядка с ошибками переходов: все, если их не больше size, иначе случайная
    выборка size различных; возвращаются только обнаружимые через контекст: список (неисправности, мутант)."""
    K = net.K
    trans = [(b, u) for b in range(K) for u in range(net.pU)]
    total = len(trans) * (len(trans) - 1) // 2 * (K - 1) ** 2
    out = []
    if total <= size:
        for (t1, t2) in itertools.combinations(trans, 2):
            for x in range(K):
                if x == net.dB[t1[0]][t1[1]]: continue
                for y in range(K):
                    if y == net.dB[t2[0]][t2[1]]: continue
                    out.append(((t1[0], t1[1], 'T', x), (t2[0], t2[1], 'T', y)))
    else:
        seen = set()
        while len(out) < size:
            t1, t2 = rng.sample(trans, 2)
            x = rng.choice([v for v in range(K) if v != net.dB[t1[0]][t1[1]]])
            y = rng.choice([v for v in range(K) if v != net.dB[t2[0]][t2[1]]])
            key = (min((t1, x), (t2, y)), max((t1, x), (t2, y)))
            if key in seen: continue
            seen.add(key)
            out.append(((t1[0], t1[1], 'T', x), (t2[0], t2[1], 'T', y)))
    res = []
    for f in out:
        M = net.mutant(f)
        if net.distinguishable(M):
            res.append((f, M))
    return res


def third_order_sample(net, rng, size, tries=20000):
    """Случайная выборка (до size) обнаружимых мутантов третьего порядка с ошибками переходов."""
    K = net.K
    trans = [(b, u) for b in range(K) for u in range(net.pU)]
    if len(trans) < 3 or K < 2: return []
    out, seen, t = [], set(), 0
    while len(out) < size and t < tries:
        t += 1
        ts = rng.sample(trans, 3)
        f = tuple(sorted((b, u, 'T', rng.choice([v for v in range(K) if v != net.dB[b][u]])) for (b, u) in ts))
        if f in seen: continue
        seen.add(f)
        M = net.mutant(list(f))
        if net.distinguishable(M): out.append(M)
    return out
