<!-- 리뷰용 한국어 번역본. 정본은 영어 파일(project_journey.md)입니다. 기술 용어·유전자명·툴명·조건명은 영어로 유지합니다. -->

# SpatialPerturb — 프로젝트 진행 여정 & 과정 기록

> 이 프로젝트가 **실제로 어떻게 흘러왔는지**를 서사로 정리한 문서입니다 — 밟아온 경로,
> 방향을 틀게 만든 문제들, 결정하기 어려웠던 지점들, 그리고 앞으로 향하는 방향.
> 실패를 제외한 `hackathon_progress.md`, 날짜별 `1-Doc/` attempt log과 짝을 이루는
> "스토리" 문서로, *무엇을* 했는지뿐 아니라 *지금 프로젝트가 왜 이런 모습인지*를
> (미래의 세션을 포함해) 누구든 이해할 수 있게 하려고 만들었습니다.
>
> EN master. Korean twin: 이 파일(`project_journey_ko.md`).
> 최종 수정: 2026-07-11.

---

## 0. 이 프로젝트가 무엇인가 (한 문단)

**SpatialPerturb**의 핵심 질문은 하나입니다: *종양 안 어디에서 T-cell perturbation이
실제로 의미 있을까?* 원래 연결하도록 설계되지 않은 두 데이터셋을 잇습니다 —
**Marson–Pritchard genome-scale CD4⁺ T cell Perturb-seq atlas** (gene knockdown이
CD4⁺ T cell에 어떤 영향을 주는지를 **Rest / Stim8hr / Stim48hr** 상태에서 측정)와
**tumor spatial transcriptomics** (T-cell state가 tumor·myeloid·stromal 이웃 대비
어디에 위치하는지). 목표는 causal claim이 아니라 **target discovery를 위한 hypothesis
generation**입니다: *target × niche* priority map을 만들되, activation context에 따라
효과의 방향이 뒤집히는 regulator가 드러나도록 **resting과 stimulated 상태를 따로
scoring**합니다. 모든 projection은 exploratory로 라벨링합니다.

---

## 1. 지금까지의 흐름 — 다섯 단계(phase)

프로젝트는 직선으로 진행되지 않았습니다. 알아볼 수 있는 다섯 phase를 거쳤고, 그중 둘은
알게 된 사실 때문에 강제된 pivot이었습니다.

### Phase 1 — 계획 & framing (Day 0, ~07-06/07)
7일 계획서(`hackathon_plan.md`, ~28 KB)를 작성했고, 특징은 **validity와 over-
interpretation risk를 맨 앞(§1)에** 놓았다는 점입니다 — 뒤에 덧붙이는 게 아니라. 계획서는
프로젝트의 추론 한계를 명시합니다: *prioritization / hypothesis-generation*으로는
valid하지만, knockout이 niche를 causal하게 remodel한다는 근거로는 **아님**. 가장 큰
과학적 리스크를 먼저 지목했고 — **cell-type confounding** ("T-cell program" spatial
신호가 T-cell *state*가 아니라 *abundance*로 생길 수 있음) — 가장 큰 방법론적 리스크도
지목했습니다 — **작은 shared gene set 위에서의 unstable transfer**.

### Phase 2 — 데이터 정찰 & 데이터셋 확정 (07-08/09)
무거운 compute 전에 "Day-1 blocker" 네 질문을 먼저 해결한 단계:
1. Perturb-seq 배포본에 condition별 precomputed signature가 있는가?
   → **Yes** (`DE_stats.h5ad`, 10,282 measured genes) — 덕분에 22M-cell raw data를
   재처리할 필요가 없었음. 1주짜리 프로젝트를 가능하게 만든 결정.
2. 어떤 spatial dataset? → **CosMx NSCLC를 primary로 확정** (§2, D1 참조).
3. Deconvolution 필요? → **아니오** — CosMx는 single-cell이며 `cell_type`이 제공됨.
4. Shared gene space가 충분히 큰가? → 그 시점 **501 shared genes**, stability floor
   위이지만 작아서 **pathway/module-level transfer를 default로** 삼고, gene-level
   score는 항상 overlap 크기를 함께 보고하도록.

### Phase 3 — 재정의(reframe): "outcome" → "어디서 의미 있는가" (07-11)
프로젝트의 중심 질문이 **사용자에 의해 재정의**되었습니다 — *perturbation의 조직 내
outcome을 예측한다*에서 *어떤 Perturb-seq target이 어떤 tumor niche에서 의미 있는지
prioritize한다*로. 단순한 표현 변화가 아닙니다. method 순서가 바뀌었습니다:
**niche definition이 Step 1이 되었고**, Moran's I / spatial-autocorrelation은
**"program을 뽑는 방법"에서 "niche가 실제로 spatially organized인지 확인하는
validation"으로 강등**되었습니다. Scoring 대상은 곱(product)이 되었습니다:
`Spatial Priority = Response Compatibility × Desired Program Change × Niche
Specificity`. 이 reframe이 지금 작업이 Moran's-I gene list가 아니라 niche를 중심으로
도는 이유입니다.

### Phase 4 — niche를 bottom-up으로 구축 (07-11, Session 1–2)
atlas에 제공된 CellCharter niche label을 그대로 믿는 대신, T-cell niche를 처음부터
만들었습니다: T cell별 **Gaussian-kernel neighborhood expression** (σ=50 µm,
3σ=150 µm cutoff, *측정된 median nearest-neighbor distance 46 µm에서 설정* — 추측한
상수가 아님) → PCA → kNN → Leiden. 결과는 20개 T-cell niche; CellCharter 대비
ARI=0.208 / NMI=0.338 — 중간 정도 일치이며, 새 niche가 reference를 **세분(subdivide)**
합니다 (TLS와 immune niche가 각각 여러 sub-niche로 쪼개짐). **묻지 않아도 밝힌 caveat**:
neighborhood expression을 batch-correct하지 않아 niche에 강한 per-sample/batch effect가
나타남.

### Phase 5 — CD4 naive로 좁히기 + state matching (07-11, Session 3)
Scope를 **의도적으로 좁혀** CD4 naive cell, LUAD-5 replicate R1–R3만으로 한정
(§2, 결정 D3 참조). 이 anchor 위에서: 21개 CD4 naive niche, 그리고 tumor CD4 naive
cell이 Perturb-seq **Rest** state에 가장 가깝다는 **state-matching** 분석 —
argmax Pearson = Rest가 cell의 **81.7%**, **21개 niche 전부**에서 (mean Pearson
Rest 0.194 > Stim48hr 0.175 > Stim8hr 0.170). 생물학적으로 예상되는 결과이며
(tumor 내 naive cell은 급성 TCR-stimulation 상태가 아님), downstream projection이
*어떤* Perturb-seq condition에 무게를 둬야 하는지를 정박시켜 줍니다.

---

## 2. 정말로 어려웠던 의사결정들

답이 자명하지 않아 (대개 사용자와 함께) 실제 고민이 필요했던 분기점들입니다.

**D1 — 어떤 spatial dataset을, 몇 개.**
요구사항 자체가 *프로젝트 중간에 바뀌었습니다*: 처음 고른 Visium HD colorectal은,
Visium HD가 bin-grid 데이터이고 **single-cell segmentation도 cell-type annotation도
없다**는 걸 알게 되자 (Space Ranger 3.0은 cell boundary 미제공) 폐기. 요구사항이
"single-cell resolution *이면서* cell-type annotation이 있을 것"으로 굳었습니다.
CosMx **breast** 옵션도 조사했으나 clean한 public download가 없어 포기. 착지점:
**CosMx NSCLC를 primary** (public, annotated, CD4 subset 존재, TLS-like 'lymphoid
structure' niche를 내장 positive control로 제공) **+ CosMx BCC를 secondary**로.
여기서의 어려움은 "명확한 immune architecture"와 "실제로 public하게 다운로드 가능하고
annotated"를 저울질하는 것 — 매력적인 여러 dataset이 두 번째 조건에서 탈락했습니다.

**D2 — niche로 분류할까, cell identity로 분류할까?** *(진짜 긴장 지점)*
핵심 개념 결정: Perturb-seq를 공간과 연결할 때 분류 단위가 **niche**(cell이 어디 있는가)
인가 **cell identity/state**(cell이 무엇인가)인가? ProjecTILs, Yeh et al. Nat Immunol
2024, SpatialProp에 근거한 결론은 **cell-identity-first**였습니다: T cell을 state로
분류하고, niche는 *spatial-context layer*(매칭된 state가 어디에 떨어지는가)로 두되
분류 단위로 쓰지 않음. 이 결정은 여전히 살아있는 조직 원리이고, Session 3가 anchor를
"TLS niche"가 아니라 "CD4 naive cell"로 잡은 이유입니다.

**D3 — scope를 얼마나 좁힐까.**
8개 sample·5개 T subtype를 다 유지할지 좁힐지. 강하게 좁혔고, 근거가 분명합니다:
**Perturb-seq는 CD4-axis 실험이며 — CD8은 애초에 profiling되지 않았습니다** — 따라서
perturbation response를 CD8 cell에 projection하는 건 정당화되지 않고, Treg는 in-vitro에서
생겨난 state로만 등장합니다. 그래서: **CD4 먼저, naive 먼저, LUAD-5 replicate만**
(cross-patient 일반화 전에 within-patient reproducibility). 비용은 분명하지만(한 histology,
한 T subtype만 봄), 이득은 downstream의 모든 주장이 Perturb-seq 데이터가 실제로 뒷받침할 수
있는 범위 안에 머문다는 점입니다.

**D4 — pathway-level vs gene-level transfer.**
CosMx panel이 whole-transcriptome Perturb-seq axis와 겨우 수백 개 gene만 공유하므로
gene-level transfer는 통계적으로 불안정합니다. 결정: **pathway/module-level transfer를
default로** (activation, cytotoxicity, exhaustion, IFN, proliferation, antigen
presentation, metabolic stress), shared gene 위의 gene-level cosine/GSEA는 항상
overlap 크기를 보고하는 *secondary* check로.

**D5 — reference-projection method: ProjecTILs vs signature scoring.**
ProjecTILs는 바로 쓸 수 있는 human CD4⁺ TIL reference atlas를 제공하지만
full-transcriptome scRNA-seq용으로 설계되어, sparse count의 960-gene imaging panel에서는
projection 품질이 떨어질 수 있습니다. **권장 main pipeline**은 CD4 identity axis에 대해
**robust signature scoring** 쪽으로 기울고, ProjecTILs를 대안으로 둡니다 — 그리고
SpatialProp/Celcomen 계열 GNN tissue-propagation은 hackathon 기간에는 너무 무거워
*의도적으로 뒤로 미룬 더 정교한 follow-up*으로 기록했습니다.

---

## 3. 문제가 우리를 문 지점들 (그리고 해결)

깔끔한 report에서는 대개 안 보이는 막다른 길과 수정 내역입니다.

**P1 — `shared genes = 0` (조용한 킬러).**
첫 state-matching 실행에서 violin plot이 **비어서(blank)** 나왔습니다. log에서 찾은
근본 원인: **CosMx `var_names`는 HGNC gene symbol** (AATK, ABL1…)인데
**pseudobulk `var_names`는 Ensembl ID** (ENSG…)이고 symbol은 `var['gene_name']`에
숨어 있었음. 단순 intersection이 비어버림. 수정: CosMx symbol을 pseudobulk
`gene_name` column을 통해 매핑 → **912 shared genes**. 노트의 헷갈리던 "501 vs 912"도
정리됨 (501은 `DE_stats`에 대한 이전 filter, 912는 pseudobulk와의 symbol intersection).
채택한 상시 규칙: **similarity를 돌리기 전에 항상 shared-gene count > 0을 검증** — 빈
gene set은 조용히 blank plot을 만들고 나중에 PCA `n_components=-1` crash를 냄.

**P2 — niche의 batch effect.**
첫 (8-sample) niche 실행에서 여러 sample이 한 niche로 붕괴 (LUAD-12 → 94.6%가 한 niche;
LUSC-6 → 88.8%가 한 niche). **batch-correct하지 않은** neighborhood expression 위에서
clustering해 sample 수준 expression 차이가 지배한 것으로 정확히 진단. 숨기지 않고 caveat로
드러냈고, LUAD-5 replicate로 좁히자(Session 3) 완화됐으나 완전히 사라지진 않음 (일부 CN
niche는 여전히 replicate-specific). Batch correction (neighborhood-expression PCA에
harmony)은 유력한 next step으로 남아 있음.

**P3 — Auditor 반려: neighbor composition vs sample distribution.**
자동 review가 niche annotation이 **neighbor cell-type composition 대신 per-sample
distribution을 대용으로** 썼다고 flag — 계획이 요구한 것과 다른 양. follow-up job에서
niche별 실제 Gaussian-weighted neighbor composition을 계산해 수정했고, 이 수정은
생물학적으로 자기 검증됨 (TLS-mapped niche가 높은 B-cell neighbor fraction을 보임 —
TLS의 B/T zone과 일치).

**P4 — 인프라 마찰.**
비과학적 문제 몇 가지가 실제 시간을 잡아먹었고 지금의 작업 방식을 만들었습니다:
- **PACE-Phoenix login node가 간헐적으로 무응답** (`echo alive`조차 timeout). 대응:
  deliverable을 artifact로 저장(서버 독립적)하고, connectivity probe가 성공할 때까지
  서버 파일 쓰기를 미룸.
- **HOME이 20 GB이고 ~77% 참** → 모든 conda env와 output을 **scratch** prefix로,
  HOME에는 절대 쓰지 않음.
- **matplotlib RGB-tuple color**를 pandas categorical `.map()`에 넘기면
  `TypeError: Expected tuple, got str` → `to_hex` + `.astype(str)`로 수정.
- **command line의 큰 base64 blob**은 "Broken pipe" (arg-length 한계) → 파일별
  heredoc으로 개별 작성.
- **한글 matplotlib label**이 AppleSDGothicNeo 등록 전에는 tofu(□)로 렌더됨.
- **번역 truncation**: LLM으로 긴 문서를 번역하면 `max_tokens`에서 문서 중간에 조용히
  잘려 bibliography 항목이 누락될 수 있음 — 이제 header가 아니라 bullet 단위로 EN-vs-KO
  item-count diff로 잡음.

---

## 4. 현재 상태 (2026-07-11 기준)

**완료:**
- plan, data-description, related-work, background literature review 문서
  (모두 bilingual EN master + KO twin).
- 데이터 다운로드 (~75 GB, scratch): Perturb-seq `DE_stats` + pseudobulk +
  supplementary table; CosMx NSCLC clustered `.h5ad`; CosMx BCC `.rds`.
- `env_spatial` 구축·검증 (scanpy 1.11.5, squidpy 1.8.2, anndata 0.12.19).
- Perturb-seq EDA / sanity check 통과 (on-target KD ~62%; TCR gene은 stim-specific;
  ~101개 Stim-specific, ~152개 Rest-specific regulator).
- CosMx atlas 시각화; T-cell-rich niche flag (lymphoid structure 3.4×,
  immune 3.1× enrichment).
- T-cell niche를 두 방식으로 구축 (8-sample all-T; CD4-naive LUAD-5).
- State matching: tumor CD4 naive ↔ Perturb-seq **Rest** 확립.

**아직 안 함 (남아 있는 진짜 과학):**
- **Step 2 본체** — spatially organized program 각각을 강화/억제하는 Perturb-seq CD4
  perturbation을 ranking, **Rest vs Stim 따로 scoring**, context-flip flagging,
  permutation null, FDR 포함.
- **Step 3** — `Spatial Priority` target × niche map 자체 (primary deliverable),
  framework figure, labeled spatial atlas.
- niche/program이 spatially organized임을 보이는 **Moran's I validation** (reframe에
  따라 이제 validation 단계).
- 선택: batch correction; BCC 재현; `pert2state_model` / CONCERT 비교; Streamlit demo.

---

## 5. 프로젝트가 향하는 방향

조직 논리는 이제 정해졌습니다: **T cell을 identity/state로 분류(CD4 axis, naive 먼저),
niche는 spatial-context layer로 사용, Perturb-seq response를 pathway level에서 transfer,
Rest vs Stim을 따로 scoring해 context-flipper를 headline output으로.** 바로 다음 수는
Step 2 — response-vector matrix를 program별·condition별 perturbation ranking으로
바꾸는 것 — 이고, 그다음이 프로젝트의 실제 deliverable인 Step 3 target × niche map입니다.
minimum viable submission(한 T subtype, 한 histology, null 포함 per-state perturbation
ranking, labeled niche × target map 하나, 정직한 limitation) 이후의 모든 것은 강화이며,
시간이 부족할 때의 cut-line은 계획서에 명시돼 있습니다.

정직한 framing은 처음부터 끝까지 유지되었습니다: 이것은 **prioritization instrument이자
hypothesis generator**이고, 모든 map은 **exploratory projection**이며, in-vitro→in-tissue
gap과 ~900-gene shared space는 limitation 문단에 묻는 게 아니라 **모든 score와 함께**
명시됩니다.
