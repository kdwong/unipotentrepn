# Orbit-to-theta-chain-to-PBP calculator

`theta2_pbp.py` combines the VALUE-path rules of the user's `theta2.py` with
Ma's BMSZ painted-bipartition implementation.  It computes paths first and
then computes their tableaux.  It does not infer, double, or otherwise alter
a painted-bipartition count.

## Run the webpage

Double-click `start_theta2_pbp_web.bat`, leave its terminal open, and use

```text
http://127.0.0.1:8000/
```

Enter the all-even orbit and click **Calculate**.  The final real form is
always `SO(n,n+1)`; the opposite orientation is not calculated.

Each `k` section shows:

1. the exact right-trivial targets of degrees `k` and `p-k`;
2. every packet-certified theta path reaching either target;
3. the raw Ma tableau attached to each path; and
4. grouping only after equality of the actual `(tau_wp,wp)` parameters has
   been checked.

Every concrete twist history names the orthogonal group at which its twist is
applied, for example

```text
O(1,2) dt -> O(4,3) tt -> O(8,9) dt
```

This compact history is shown directly below `Path X` in the webpage summary;
the full lowest-harmonic chain remains available when the path is expanded.
Thus the last entry of a history ending in `O(8,9) dt` appears in that chain as
`final twist[dt]`; a grouped list such as `final twist[dt,tt]` is not used for
a concrete realization.  Orthogonal lowest harmonics are written in padded
type notation, for example

```text
O(6,7)(111000|0000000) [untwisted lowest harmonic]
final twist[dt] -> O(6,7)(111000|0000000)
```

rather than as the internal VALUE lists `L=[1,1,1] R=[]`.

For example, for `O^vee=(2,2,2,2,2,2)`, final form `SO(6,7)`, and `k=1`, the
two exact endpoint types are

```text
(100000|0000000)
(111110|0000000)
```

Ten certified paths reach these two endpoint types and are grouped into two
different PBPs.  The page obtains both from Ma's calculation; it does not
manufacture four representations by multiplying counts.

## Run the Python report

Edit the one orbit line near the beginning of `theta2_pbp.py`:

```python
DUAL_ORBIT = (6, 4, 2, 2, 2)
```

Then run:

```powershell
Set-Location 'C:\Users\kayue\Documents\BMSZ'
.\.venv\Scripts\python.exe theta2_pbp.py
```

The final form is fixed to `O(n,n+1)`.  To override only the orbit from the
command line, use for example:

```powershell
.\.venv\Scripts\python.exe theta2_pbp.py --orbit 6,4,2,2
```

For this example the program therefore uses `O(7,8)`.  In general,

\[
 (p,q)=(n,n+1),\qquad 2n=|O^\vee|.
\]

Add `--show-diagrams` to print the tableaux vertically.

Running `python theta2_pbp.py` with another interpreter also works: the file
relaunches itself inside the prepared `.venv`, which contains the required
`multiset` package.

## What a `k` bin means

Fix the ordered final form `O(p,q)`.  The requested right-trivial exact types
are

\[
 \bigwedge^k\mathbb C^p\boxtimes\mathbf 1,
 \qquad
 \bigwedge^{p-k}\mathbb C^p\boxtimes\mathbf 1.
\]

The program constructs both.  Only after their theta paths and PBPs are known
does it put them under the common connected highest-weight label `1^k`.
When `p=2k`, the two exact exterior degrees coincide, but every concrete
character-twist history is still expanded before PBP deduplication.

The final endpoint printed on each path is full-dimensional.  Thus a path in
the `k=1` group for `O(6,7)` ends in either
`(100000|0000000)` or `(111110|0000000)`, rather than hiding both behind the
rank-compressed label `(100|000)`.

The right-determinant representatives are not separately listed because the
requested endpoint has right side zero.  Multiplication by the total
determinant pairs them with the displayed right-trivial representatives and
does not create a new SO painted bipartition.  The internal BMSZ outer
`epsilon` is retained as audit metadata, but it is not turned into a second
PBP count.

## Calibration: `O^vee=(2n)`

For a one-row orbit the tower is the direct lift

\[
 Mp(0)\longrightarrow O(2n+1).
\]

For every `1 <= n <= 50`, the regression suite compares the program with the
complete Ma final packet for `O(n,n+1)`.  In all fifty cases:

- only `k=0` is nonempty;
- the two exact degrees are `0` and `p`;
- the twists are respectively `tt` and `dt`;
- the two returned raw tableaux equal the complete set of Ma PBPs of that
  real form; and
- every higher-`k` bin is empty.

For the smallest nontrivial example, `O^vee=(4)` and `SO(2,3)`, this gives
exactly

\[
 \varnothing\times\begin{matrix}s\\d\end{matrix}\times B^+ ,
 \qquad
 \varnothing\times\begin{matrix}s\\r\end{matrix}\times B^- .
\]

They come from the explicit final twists `tt` and `dt`; the second is not an
arithmetically supplied determinant copy of the first.

## Raw `PBP(O^vee,wp)` shapes

For `wp=empty` and `O^vee=(r_1,r_2,...)`, the reference type-B column shape is

\[
 P=(r_2/2,r_4/2,\ldots),\qquad
 Q=(r_1/2,r_3/2,\ldots).
\]

A nonempty `wp` changes this shape.  The program retains the domain element
`tau_wp` in

\[
 \bigsqcup_{\wp}\operatorname{PBP}(O^\vee,\wp).
\]

Ma's bijection with the `wp=empty` reference set is used only to recover and
certify `wp`; it is not used to replace the displayed tableau.  Thus for
`O^vee=(6,4,2,2,2)` the program keeps

```text
wp=empty:    (2,1 | 3,1,1)
wp={(2,3)}: (1,1 | 3,2,1)
```

and the second entry remains the user's `(2 x 321)` example.

## Derived theta tower

The program derives all group dimensions from the orbit.  It forms reverse
cumulative row sums, appending a zero row first when the number of rows is
even.  A type-B stage of suffix size `S` is `O(S+1)`, while a type-M stage is
`Mp(S)`.

For `O^vee=(6,4,2,2)` this gives

```text
theta2 totals: (1,2,5,8,15)
Mp(0) -> O(1) -> Mp(2) -> O(5) -> Mp(8) -> O(15).
```

Under the fixed convention, the last stage is `O(7,8)`.  The program never
also calculates `O(8,7)`.

## The counting rule is checked, not imposed

Let `m_k` be the number of row-submultisets of `O^vee` whose rows sum to `2k`.
The proposed rule is

\[
 N_k=\begin{cases}
 m_k,&n=2k,\\
 2m_k,&n\ne2k.
 \end{cases}
\]

This rule is compared with the output only after the calculation.  It is not
used to create paths or PBPs.  For the fixed `SO(n,n+1)` convention,
the rule agrees in all `261/261` fine-K bins for the 66 all-even orbits through
total orbit size `16`.  This is a bounded computational check, not a proof of
the rule in arbitrary rank.  Every displayed number still comes solely from
certified paths and deduplication of their raw Ma PBPs.

## Verification commands

Focused examples, the complete one-row calibration, and raw `wp` shapes:

```powershell
.\.venv\Scripts\python.exe tests\test_theta2_pbp.py
```

Direct comparison with the user's original `theta2.py` for all 58 multi-row
all-even orbits through total size `16`:

```powershell
.\.venv\Scripts\python.exe tests\test_theta2_reference_parity.py
```

Web API and real local HTTP requests:

```powershell
.\.venv\Scripts\python.exe tests\test_theta2_pbp_web.py
```

Structural replay of every intermediate/final twist and every Ma PBP shape
for all 66 all-even orbits through total size `16`:

```powershell
.\.venv\Scripts\python.exe tests\test_theta2_pbp_exhaustive.py --max-total 16
```

Large regression for `O^vee=(8,6,4,4,2,2,2,2)`:

```powershell
.\.venv\Scripts\python.exe tests\test_theta2_pbp_86442222.py
```

The fixed five-chain checker for the originally supplied `O(7,8)` examples
remains available separately:

```powershell
.\.venv\Scripts\python.exe tests\test_theta_chains.py
```
