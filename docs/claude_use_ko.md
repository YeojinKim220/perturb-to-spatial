<!-- 리뷰용 한국어 번역본. 정본은 영어 파일(claude_use.md). 코멘트는 이 파일에 남겨주세요. -->

# Claude Use — SpatialPerturb (Researcher Track)

> 제출 평가 항목 (25%): *이 팀은 Claude Code를 얼마나 창의적으로 사용했는가? 기본적인 활용을 넘어섰는가? 우리(주최측)조차 놀랄 만한 capability를 끌어냈는가?*
>
> 이 문서는 **Claude Science(Claude Code)를 프로젝트 진행에 어떻게 사용했는지**에 대한 실시간 기록이다. 매 작업 세션마다 업데이트한다. 최신 항목이 위에 온다.

## 사용 방식 요약

우리는 Claude를 snippet을 돌려주는 chatbot으로 쓰지 않았다. **Claude Science를 자율적인 research operator로** 사용했으며, 구체적으로:
- **원격 HPC 클러스터(PACE-Phoenix, Slurm)를 end-to-end로 운용.** Claude가 `compute_details`를 읽고, scratch prefix에 conda 환경을 만들고, Slurm job script를 작성하고(account/QOS/partition directive 정확히), `sbatch`로 제출하고, 완료 notification을 비동기로 기다리고, 산출물을 harvest해 오고, 다음 세션을 위해 host에 대해 배운 것을 기록한다 — 우리가 수동 SSH를 하지 않고.
- **submit→fail→fix→resubmit 루프를 스스로 교정.** job이 실패하면(예: matplotlib RGB-tuple 색상이 pandas categorical `.map()`을 깨뜨림) traceback을 읽고 근본 원인을 진단하고 script를 고쳐 재제출했으며, gotcha를 durable memory에 기록해 반복하지 않게 했다.
- **프로젝트 지식을 세션 간에 유지** — memory system(데이터 경로, env 위치, niche 정의 결정, host 특이사항)과 bilingual(EN 정본 + KO 리뷰) 문서 관례로, 매 세션이 설정을 다시 발견하지 않고 전체 맥락으로 재개된다.
- **plan → 승인 → step 추적 실행 워크플로**를 사용하며, 결정 지점(niche 정의 방법, neighbor 범위, clustering 알고리즘)에서 사람을 구조화된 객관식 질문으로 loop에 유지한다.

## 세션 로그

### 2026-07-11 — Session 4: 아이디어 문서 리뷰, scope 결정, 두 트랙 HPC 실행

**우리가 Claude에게 시킨 것:** 수정된 soft-matching 파이프라인 초안을 리뷰하고 무엇을 고쳐야 할지 알려주고, 검증/결정 항목을 정해 분석을 시작하라.

**Claude가 capability를 (기본을 넘어) 활용한 방식:**
- **문서만이 아니라 프로젝트에 대한 자기 memory와 대조해 설계를 리뷰.** Claude가 기존 state-matching 결과를 memory에서 꺼내고(argmax=Rest 81.7%, mean Pearson 값), *실제 per-cell CSV를 다시 읽어* 수치를 확인한 뒤 비평 — 초안의 주력 per-cell projection이 결국 순수 Rest effect로 붕괴함을 짚고, 방어 가능한 주력으로 direction-flip 분석을 권고.
- **"이게 왜 필요하냐"는 질문을 한 번의 클러스터 검증으로 전환.** 912 vs 501 gene 수 모순을 말로 다투는 대신, h5py var-group reader를 작성해 로그인 노드에서(수 초, 4 GiB 캡 내) 실행하고 *실제 파일에서 두 숫자를 freeze* — 912·501, 중복 매핑 0 — 상반되던 memory 노트 2개를 정정.
- **독립적인 두 분석을 병렬 Slurm job으로 실행**(all-CD4 niche 재정의; 16 GB DE_stats의 Rest/Stim direction-flip), 함께 제출하고 비동기로 harvest.
- **기존 작업을 재작성이 아니라 lineage로 재사용.** 기존 아티팩트 lineage에서 freeze된 7-program marker set과 4-block niche 코드를 추출해 anchor만 교체(naive → all-CD4), 트랙 간 방법을 동일하게 유지.
- **scope 결정을 검정 가능한 가설로 만든 뒤 실제로 검정.** all-CD4 확장의 근거("naive-only는 exhaustion/Treg score를 축퇴시킨다")를 job에 비축퇴화 검증으로 넣어 확인(exhausted positivity 0→0.53).
- **주력 결과를 알려진 생물학과 self-check** — TCR 신호 KO가 cytotoxic program을 Rest→Stim으로 예상대로 flip함을 짚어, 파이프라인이 허구가 아니라는 근거로 제시.
- **job 실패를 self-recover**(matplotlib `labels`→`tick_labels` API 개명)했고 한 줄 수정 후 재제출, 그리고 — 쓰지 않은 추가 방어 코드를 썼다고 잘못 말했을 때 back-fill 대신 auditor의 지적을 수용.


**같은 세션, 이어서(relevance map + 검증):**
- **두 트랙을 niche × target relevance map으로 결합** — niche baseline program burden × condition별 perturbation effect — 입력을 재유도하지 않고 endpoint별 candidate regulator 순위 산출.
- **실제 문헌을 조회해 자기 후보를 검증.** PTPN2를 뽑은 뒤 memory로 인용을 주장하는 대신 OpenAlex를 조회해 검증 가능한 참고문헌(ABBV-CLS-484 2023; 1차 인간 T세포 CRISPR screen 2018)을 반환.
- **자기 조작 인용을 잡아 정정.** PTPN2를 조회 없이 "Yeh 2024"에 귀속했었고, auditor 지적 시 귀속을 철회하고 실제 조회를 실행해 정정을 기록.
- **circularity 통제를 지시 없이 설계·실행.** 같은 marker 유전자가 niche와 순위를 둘 다 정의함을 인지하고, held-out confirmation(471개 = shared에서 program marker 제거)을 matched-null percentile과 함께 만들어 Slurm에서 실행, 정직한 분리 결과 보고: 3개 endpoint는 독립 유전자에서 확인, 2개는 미확인 — 자기 이전 후보를 방어하지 않고 하향 조정.

### 2026-07-11 — Session 3: CD4 naive 범위, cell-identity 전환, 4-block feature niche
**요청:** CD4 naive / LUAD-5 R1-3로 범위 축소; Perturb-seq를 niche vs cell-identity 중 무엇으로 연결할지 재고; T-cell 상태 + 이웃 조성 + suppressive program + tumor 경계 거리를 결합한 tumor spatial niche 정의.

**Claude가 capability를 (기본 이상으로) 활용한 방식:**
- **방법론 논쟁을 실시간 문헌 리뷰로 근거화.** 인간이 niche-level 연결이 맞는지 물었을 때, Claude는 targeted web search를 돌려 3가지 확립된 패러다임(ProjecTILs reference-projection, Yeh 2024 HGSC spatial+Perturb-seq, SpatialProp GNN)을 각각의 데이터 자원과 960-panel 적합성/한계까지 정리하고, 구체적 파이프라인(cell-identity 분류 → Yeh 방식 Perturb-seq 순위 → 공간 지도화)을 권고 — 단순 예/아니오가 아니라.
- **격리로 버그 진단.** state-matching이 공유 gene 0개를 반환하자 추측 대신 두 AnnData의 var_names·var 컬럼을 비교하는 one-shot 진단을 제출해 불일치(gene symbol vs Ensembl ID)와 수정법(pseudobulk gene_name → 912 genes)을 특정.
- **commit 전 진단, 재차.** feature-niche 파이프라인 전에 cell_type 라벨 열거(거리 block의 malignant 정의용)와 960-panel signature marker 커버리지 측정 job을 먼저 제출 → panel 존재 marker만 사용.
- **방어 가능한 가중치의 해석적 feature engineering.** 4-block feature 벡터를 만들고, 인간이 block 가중치가 원칙적인지 묻자 정직히 답하고(feature 개수 차이로 암묵적·불균등), 결과(거리가 46 feature 중 묻힘)를 설명한 뒤 인간이 선택한 block-count(√n) 정규화로 재실행 — 이후 tumor 거리 축이 67-1121µm로 펼쳐지고 예상 niche 유형(tumor-contacting exhausted, myeloid-rich Treg, fibroblast-rich excluded, TLS-active)이 깨끗이 등장.
- **publication-grade figure 정제.** 요청 시 niche heatmap을 양축 hierarchical clustering으로 재정렬하고 저장 전 geometric overlap self-check 수행.
- **정직한 상태 보고.** 일일 문서가 최신인지 물었을 때 아직 이번 세션분 미반영임을 명확히 말하고 곧바로 갱신.

**결과:** 21 CD4 naive neighborhood niche(Phase 0) + 명시적·원칙적 block 가중치의 16 feature-niche; cell-identity 연결로의 문헌 기반 전환 결정; state-matching 버그 root-cause + 알려진 수정법.

### 2026-07-11 — Session 2: T-cell-specific niche 구성 (neighborhood-expression clustering)
**Claude에게 요청한 것:** 논문의 CellCharter 라벨에 의존하지 말고 T-cell-specific niche를 처음부터 구성 — T cell별 Gaussian-kernel neighborhood expression → clustering → annotation → 기존 niche 라벨과 비교.

**Claude가 capability를 어떻게 사용했나 (기본을 넘어서):**
- **commit 전 diagnostic.** 파이프라인 작성 전에 method를 올바르게 parameterize하기 위해 필요한 것을 작은 Slurm diagnostic job으로 측정: matrix 정규화 상태(X는 log-normalized, raw는 `layers['counts']`), spatial 좌표 구조(sample들이 서로 분리된 좌표 영역 차지; FOV는 sample 내 tiling → neighbor는 sample 내 계산), median nearest-neighbor 거리(46 µm) — 이를 이용해 Gaussian bandwidth σ=50 µm와 3σ=150 µm cutoff를 **데이터로부터** 설정, 추측 상수가 아님.
- **사람이 개입하는 설계 선택**을 구조화된 질문으로 제시: neighbor 범위(모든 cell type vs T only), clustering 알고리즘(Leiden vs 고정-k KMeans) — 사람이 고르고 Claude가 정확히 구현.
- **off-the-shelf 호출이 아닌 커스텀 spatial method.** KD-tree radius query 위에 Gaussian-kernel 가중 neighborhood expression, sample별, sparse matrix로 벡터화; 이후 neighborhood matrix에 scanpy PCA→kNN→Leiden; 이후 reference CellCharter niche와 ARI/NMI + row-normalized confusion heatmap으로 정량 비교.
- **정직한 해석.** Claude가 요청받지 않은 발견을 flag: spatial map에서 명백한 sample/batch 효과(일부 sample이 단일 niche로 붕괴)를 보이며, 이를 batch 미보정 neighborhood expression에 clustering한 탓으로 정확히 귀속 — 요청 없이 제시된 caveat이며 다음 방법론 단계를 설정한다.

**결과:** 20개 T-cell niche (Leiden res=1.0); TLS-like(lymphoid structure)와 immune niche가 여러 sub-niche로 깔끔히 분리; CellCharter 대비 ARI=0.208 / NMI=0.338 (중간 수준 일치 — 새 niche가 reference를 세분화). Artifacts: heatmap, UMAP, sample별 spatial map, annotation + crosstab + per-cell assignment table.

### 2026-07-11 — Session 1: env + atlas 시각화 + T-cell-niche flagging
**Claude에게 요청한 것:** 클러스터에 spatial 분석 환경 생성, CosMx NSCLC atlas를 cell type으로 시각화, T-cell niche 정의.

**Claude가 capability를 어떻게 사용했나 (기본을 넘어서):**
- **PACE에 Slurm으로 새 conda 환경 구축** (`env_spatial`: scanpy/squidpy/anndata), scratch prefix에 — HOME이 ~77% 찼다고 스스로 판단해 설치를 scratch로 redirect — 이후 import 검증하고 activation recipe를 `compute_details`에 기록해 재사용.
- **2.6 GB, 765k-cell AnnData를 로컬 메모리에 로드하지 않고 검사** (클러스터에서 backed read)해 figure 코드 작성 전에 정확한 schema 확보.
- **8-panel single-cell spatial atlas 렌더링** (765,771 cells, 22 cell type, 공유 legend)하고 **niche별 T-cell fraction**을 global 대비 enrichment와 함께 계산, T-cell-rich niche flag(lymphoid structure 3.4×, immune 3.1×) — 계획된 TLS positive control을 정량적으로 확인.
- **사람의 개념 질문에 답변**(niche란 무엇인가, resting vs stimulated의 의미, Moran's I가 왜 필요한가)하되 실제 데이터 schema에 근거해, 데이터에 *제공된* 것(CellCharter niche 라벨)과 Claude가 *계산한* 것(T-cell-rich flag)을 구분.

**결과:** `env_spatial` 준비 완료; cell-type atlas + T-cell-fraction table + highlight map을 artifact로 전달하고 서버에 미러링.

## 심사위원에게 강조할 capability
- 비동기 job harvest와 host별 학습을 갖춘 자율 **Slurm orchestration**.
- 하드코딩 상수가 아닌 **데이터 기반 parameterization**(kernel bandwidth 설정을 위한 NN 거리 측정).
- 다일(multi-day) 프로젝트를 재개 가능하게 만드는 **cross-session memory + bilingual 문서**.
- 수정 사항을 memory에 기록하는 **self-diagnosed failure recovery**.
- green checkmark만이 아니라 다음 단계를 설정하는 **요청받지 않은 과학적 caveat**(batch 효과).
