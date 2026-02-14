# 最大流

N頂点M辺の最大流を、Dinic法でO(N^2 M)で求めます。    
以下、辺容量の型をTとします(使用上はcodonが自動で型推論します)。   
**すべての関数を通して、Tにあたる型は統一してください。**

## 更新履歴

- 2026/02/15
**`min_cut`に致命的なバグがあったため修正しました。**  
具体的には残余グラフではなく、現在のグラフ上での到達判定を行っていました。  
大変申し訳ありませんでした。  

## 使い方

`MFGraph(N: int) -> None`
- N頂点0辺の有向グラフを作成します。
- 制約: 0 ≤ N < $2^{30}$

`add_edge(now: int, nxt: int, cap: T) -> None`
- 頂点nowから頂点nxtに、容量capの辺を追加します。
- 制約: 0 ≤ cap

`get_edge(i: int) -> tuple[int, int, T, T]`
- **0-indexedで**i番目に追加した順辺の辺情報を取得します。
- `(now: int, nxt: int, cap: T, flow: T)` のタプルで返します。  
これはi番目の辺が頂点nowから頂点nxtに向かう容量capの辺で、now → nxtの方向に流量flowが流れていることを示します。  
ここで、0 ≤ now, nxt < N, 0 ≤ flow ≤ cap を満たします。
- 制約: 0 ≤ i < M

`edges() -> list[tuple[int, int, T, T]]`
- これまでに追加した順辺の本数をMとして、i = 0, 1, ･･･, M - 1 の順に、**0-indexedで**i番目に追加した順辺の辺情報を列挙します。
- i番目のリストの要素は `(now: int, nxt: int, cap: T, flow: T)` のタプルです。    
これはi番目の辺が頂点nowから頂点nxtに向かう容量capの辺で、now → nxtの方向に流量flowが流れていることを示します。  
ここで、0 ≤ now, nxt < N, 0 ≤ flow ≤ cap を満たします。

`change_edge(i: int, new_cap: T, new_flow: T) -> None`
- **0-indexedで**i番目に追加した順辺の容量をnew_capに、流量をnew_flowに変更します。
- 制約: 0 ≤ i < M, 0 ≤ new_flow ≤ new_cap

`flow(St: int, Gl: int) -> T`  
`flow(St: int, Gl: int, flow_limit: T) -> T`  
`flow(St: int, Gl: int, flow_limit: T, permissible_error: T) -> T`
- 頂点Stから頂点Glにフローを流し、流した量を返します。
- `flow(St, Gl)`では、頂点Stから頂点Glにフローを**流せるだけ**流します。
- `flow(St, Gl, flow_limit)`では、フローを**flow_limitを流量上限として**流します。    
追加制約: 0 ≤ flow_limit
- `flow(St, Gl, flow_limit, permissible_error)`では、**残容量がpermissible_error以下の辺を残容量0とみなし**、flow_limitを流量上限としてフローを流します。    
追加制約: 0 ≤ flow_limit, 0 ≤ permissible_error
- 制約: 0 ≤ St,Gl < N, St ≠ Gl, 答えは型Tの表現範囲に収まる
- 計算量: O(N^2 M)

`flow_capacity_scaling(St: int, Gl: int, flow_limit: T) -> T`
- **辺容量の型Tが整数型の場合に限りますが**、容量スケーリングに対応しています。
- 頂点Stから頂点Glに、flow_limitを流量上限としてフローを流し、流せた量を返します。
- 制約: 0 ≤ St,Gl < N, St ≠ Gl, 答えは型Tの表現範囲に収まる, 0 ≤ flow_limit, Tは(`int`, `Int[N]`, `UInt[N]`)のどれか
- 計算量: 辺の最大容量をCとして、O(NMlogC)

`min_cut(St: int) -> list[bool]`  
`min_cut(St: int, permissible_error: T) -> list[bool]`
- 残余グラフ上で、頂点Stから各頂点に到達可能か判定します。
- `min_cut(St, permissible_error)`では、残容量がpermissible_error以下の辺を残容量0とみなして判定します。  
追加制約: 0 ≤ permissible_error
- 計算量: O(N + M)
