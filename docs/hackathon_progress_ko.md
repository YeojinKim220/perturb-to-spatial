<!-- 리뷰용 한국어 번역본. 정본은 영어 파일(hackathon_progress.md). 코멘트는 이 파일에 남겨주세요. -->

# Hackathon Progress — 리포트 초안

**진행/결과/방법만** 담는 실시간 요약(실패·막다른 길은 제외 — 그건 `1-Doc/`의
날짜별 로그에 있음). 최종 리포트의 작업 초안이며 계속 업데이트한다.

## 목표
T세포 내재적 regulator가 공간적 면역 program에 **activation-context 특이적으로**
작용하는지, 그리고 공간 구조가 각 regulator가 작용할 **위치**를 예측하는지 검정.
target discovery를 위한 가설 생성 — 인과 주장은 아님.

## 데이터
- **Perturbation:** Marson–Pritchard genome-scale CD4+ T cell Perturb-seq
  atlas(CRISPRi, donor 4명, 조건 Rest/Stim8hr/Stim48hr). 저자 제공
  `GWCD4i.DE_stats.h5ad`(perturbation × gene logFC/zscore),
  `GWCD4i.pseudobulk_merged.h5ad`(NTC + targeting pseudobulk).
- **Spatial:** CosMx NSCLC(He 2022), 765,771 cells × 960-gene panel, 8 sample.
  주 substrate: LUAD-5 R1-3의 **all CD4 T cell**(naive + memory + Treg).

## 설계 결정
- **Cell-identity 우선, niche는 context.** T세포는 내재적 state(7 program)로
  분류하고, 공간 niche는 매칭된 state가 **어디**에 있는지를 말하는 layer일 뿐
  분류 단위가 아니다.
- **Scope = all CD4**(naive-only 아님): naive만으로는 exhaustion/Treg/Th1
  endpoint가 축퇴(score ≈ 0)되어 substrate를 확장.
- **주력 결과 = Rest vs Stim direction-flip.** per-cell soft context projection은
  exploratory 전용 — per-cell state-matching이 Rest로 붕괴하기 때문(argmax=Rest
  81.7%).
- **두 개의 gene overlap 분리 유지:** state-matching은 912개
  (CosMx ∩ pseudobulk), gene-level/held-out 확인은 501개
  (CosMx ∩ DE_stats, ⊂ 912). 모든 score에 overlap 크기 병기.

## 방법
1. **Gene-ID 정합.** CosMx HGNC symbol을 Perturb-seq `var['gene_name']`에 매핑;
   912(state-matching)·501(confirmation) shared gene freeze. 중복 symbol→Ensembl
   매핑 0개.
2. **공간 niche(all-CD4).** 세포별 4-block feature vector — (1) 7개 내재적 state
   signature(`score_genes`), (2) Gaussian-weighted 이웃 cell-type 구성(150µm),
   (3) 이웃 suppressive/ligand program, (4) 최근접 악성세포까지의 log 거리 — 각
   block z-score 후 sqrt(n) 정규화로 동등 기여; PCA → kNN → Leiden.
3. **Perturbation program effect.** DE_stats logFC에서 E_{g,c,program} =
   해당 program marker들의 mean logFC(target 유전자 제외), 조건
   c ∈ {Rest, Stim8hr, Stim48hr}, QC 통과 perturbation
   (on-target 유의, distal off-target 없음, target 발현 충분)에 대해.
4. **Direction-flip 분류.** (target, program)별로 E_Rest vs E_Stim48hr 비교:
   direction_flip(부호 반전, 양쪽 |E|≥0.10), rest_specific, stim_specific, stable.

## 결과
- **Gene overlap freeze:** 912(state-matching), 501(gene-level, ⊂ 912).
- **All-CD4 substrate 검증(비축퇴화):** 25,493개 CD4 세포에서 program score가
  실재 — fraction-positive exhausted 0.53, cytotoxic 0.35, Treg 0.35, Th1 0.27
  (naive-only에선 ≈0). Treg 세포가 Treg(1.60)·exhausted(0.59) score 최고.
- **all-CD4 feature-niche 19개.** 성숙 state가 독자 niche로 분리:
  FN1 = Treg(86.8%), FN2 = memory(80.2%); niche는 program축(Tfh-high FN10;
  Th1-high FN4/FN15)·tumor 경계 거리로도 분화.
- **Direction-flip 스크리닝(주력):** perturbation 7,828개 QC 통과;
  target×program 54,796쌍 중 **5,125개가 Rest↔Stim48hr direction-flip**,
  Tfh(1266)·cytotoxic(980)·Treg(817)·exhausted(719)에 집중.
  Positive control: TCR 신호 KO(ZAP70, CD3E, VAV1, PTPRC)가 cytotoxic program을
  Rest(음)→Stim48hr(양)으로 flip — 알려진 활성화 생물학과 일치.

- **Niche × target relevance map(트랙 결합):** Relevance(niche, target,
  endpoint) = max(0, niche program burden) × max(0, direction·E_Rest). endpoint별
  대상 niche = 부담 최대 niche(exhaustion FN10, Treg FN1, Th1 FN4, activation
  FN10, cytotoxic FN4). positive relevance 586,476건.
- **Held-out gene confirmation(circularity 통제):** niche↔perturbation cosine을
  471개 유전자(501-shared에서 program marker 30개 제거)로 재계산, QC 통과 Rest
  perturbation 6,734개 matched null과 비교 → activation(null-pct 0.93),
  Th1(0.89; PTPN2 0.98), cytotoxic(0.86) endpoint는 독립 유전자에서 **확인**,
  exhaustion(0.14)·Treg(0.39)는 **확인 안 됨**(marker 구동 → 하향 조정).
- **문헌 근거(OpenAlex 조회):** PTPN2(ABBV-CLS-484 2023; 1차 인간 T세포 CRISPR
  screen 2018), BHLHE40(CD8 exhaustion CRISPR 2023). 가장 방어 가능한 단일 hit =
  **PTPN2**(문헌 + held-out).

## 상태
- 2026-07-08: 문서 구조 확립.
- 2026-07-11: niche 파이프라인 구축; scope를 all-CD4로 설정; 912/501 freeze;
  Track A(niche + 검증)·Track B(direction-flip, 주력) 완료.
- 2026-07-11(계속): 트랙 결합 relevance map; held-out 확인 실행
  (Th1/activation/cytotoxic 유지, exhaustion/Treg 미확인); top candidate 문헌
  확인(PTPN2 최강).
- **다음:** 확인된 후보(PTPN2 등)의 direction-flip 궤적 그림; exhaustion/Treg
  endpoint를 marker-free niche 정의로 재접근; cross-patient 검증 미해결
  (단일 환자 LUAD-5).
