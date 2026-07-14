<!-- 리뷰용 한국어 번역본. 정본은 영어 파일(background_literature_review.md). 코멘트는 이 파일에 남겨주세요. -->

# CD4⁺ T 세포 Perturb-seq을 종양 공간 전사체학(spatial transcriptomics)에 연결하기 — 배경 문헌 리뷰

*"Build From the Bench" 프로젝트를 위해 작성됨: 게놈 규모 CD4⁺ T 세포 Perturb-seq atlas(Marson–Pritchard 2025)를 종양 공간 전사체학(CosMx NSCLC)과 연결하여, 활성화 맥락(activation-context)에 특이적인 방식으로 T 세포 내재적 조절인자가 작용할 위치의 우선순위를 매기고자 한다. 이는 가설 생성 도구이지 인과적 또는 예측적 주장이 아니다 — 아래의 문헌은 이 파이프라인의 각 단계가 무엇에 근거하고 있으며, 해당 분야의 방법론이 여전히 논쟁 중인 지점이 어디인지를 보여주도록 구성되었다.*

---

## 범위와 읽는 방법

이 프로젝트는 서로 대체로 분리된 세 가지 문헌에서 나온 세 가지 작업을 연쇄적으로 연결한다: 풀링(pooled) 단일세포 CRISPR 스크린(**perturbation → response-vector** 자원), 조직에 대한 공간 자기상관(spatial-autocorrelation) 분석(**프로그램이 어디에 위치하는가**), 그리고 이 둘을 용접하는 시그니처 정렬/투영(projection) 단계이다. 이 리뷰는 검증된 참고문헌 집합(2005–2026년에 걸친 77편의 1차 논문, PubMed 및 OpenAlex 대조로 DOI 확인)을 바탕으로 작성되었으며, 시간순이 아니라 이러한 방법론적 축을 중심으로 구성되었다. 그 이유는 프로젝트가 마주하는 설계 결정들 — 어떤 perturbation 자원을 쓸지, 어떤 공간 변이성(spatial-variability) 통계량을 쓸지, 어떤 정렬(alignment) 지표를 쓸지, 아티팩트를 어떻게 방어할지 — 각각이 현재 활발히 논쟁 중인 특정 방법론적 쟁점에 대응하기 때문이다. 어떤 주장이 프리프린트(preprint)에 근거하는 경우에는 그렇게 표시하였으며, 방법론이 논쟁 중인 경우에는 양측의 입장을 모두 제시하였다.

---

## 1. Perturbation 자원: Perturb-seq에서 맥락-해상(context-resolved) atlas로

풀링 단일세포 CRISPR 스크리닝은 CRISPR 라이브러리를 단일세포 전사체 readout과 결합하여, 각 세포가 어떤 유전자가 교란(perturb)되었는지와 그 게놈 전반의 결과를 동시에 보고할 수 있도록 한 2016–2017년의 네 편의 논문에 의해 거의 동시에 확립되었다. Perturb-seq 자체는 확장 가능한 바코드 설계와, perturbation을 고차원 발현 반응으로 보는 분석적 프레이밍을 도입하였다([Dixit 2016](https://doi.org/10.1016/j.cell.2016.11.038)); ([Adamson 2016](https://doi.org/10.1016/j.cell.2016.11.048))의 다중화(multiplexed) CRISPR 스크리닝 플랫폼과 ([Datlinger 2017](https://doi.org/10.1038/nmeth.4177))의 CROP-seq guide-capture 전략은, CRISPR-Cas9이 인간 세포에서 풀링 손실기능(loss-of-function) 스크린을 구동할 수 있음을 보인 이전의 시연([Wang 2013](https://doi.org/10.1126/science.1246981))을 토대로, 게놈 규모 스크린을 가능하게 만든 guide-assignment 및 라이브러리 복잡도 문제를 해결하였다. 이 프로젝트가 입력 데이터를 소비하는 방식과 직접 관련된, 반복적으로 등장하는 방법론적 흐름은 다음과 같다: 원본 단일세포 CRISPR 데이터는 노이즈가 많으며, 세포당 카운트를 신뢰할 수 있는 perturbation당 효과로 전환하려면 명시적인 통계 모델이 필요하다는 점이다. ([Duan 2019](https://doi.org/10.1038/s41467-019-10216-x))는 이를 perturbation 효과에 대한 모델 기반 추정으로 정식화하였고, atlas가 사용하는 pseudobulk-DESeq2 파이프라인은 이 원칙의 현재 표준적 구현이다.

이후 이 분야는 이 프로젝트가 물려받은 두 가지 축을 따라 확장되었다. 첫째는 규모 자체와 readout의 풍부함이다: direct guide-capture는 조합적(combinatorial) 스크린을 가능하게 하였고([Replogle 2020](https://doi.org/10.1038/s41587-020-0470-y)), 게놈 규모 Perturb-seq은 수만 개의 perturbation에 걸쳐 정보가 풍부한 genotype–phenotype 지형을 매핑하였다([Replogle 2022](https://doi.org/10.1016/j.cell.2022.05.013)) — 이는 "perturbation × 유전자" 효과 행렬 자체를 분석 가능한 대상으로 다루는 방법론적 표준 틀이다. 이러한 perturbation atlas 전반에 걸친 전사체 규모의 차등발현(differential-expression) 모델링은 여전히 활발한 연구 주제이다([Nadig 2025](https://doi.org/10.1038/s41588-025-02169-3)). 둘째, 이 프로젝트에 더 중요한 지점으로, 스크린이 1차 인간 T 세포와 생리학적으로 의미 있는 맥락(context)으로 이동했다는 점이다. 1차 인간 T 세포에서의 게놈 전체 CRISPR 스크린은 핵심 증식 및 사이토카인 조절인자가 오직 올바른 자극(stimulation) 맥락에서만 발견 가능함을 처음으로 밝혔다([Shifrut 2018](https://doi.org/10.1016/j.cell.2018.10.024)); 소진(exhaustion) 중심 스크린은 기능장애 상태의 크로마틴 조절인자를 식별하였고([Belk 2022](https://doi.org/10.1016/j.ccell.2022.06.001)); 후성유전적 편집(epigenetic-editing) 및 크로마틴 리모델러 스크린은 단일 유전자가 아니라 조절 계층구조를 규명하기 시작하였다([Otto 2023](https://doi.org/10.1016/j.molcel.2023.03.013); [Pacalin 2024](https://doi.org/10.1038/s41587-024-02213-3)).

기준(anchor) 데이터셋은 이러한 궤적의 종착점에 위치한다. Marson–Pritchard 게놈 규모 CD4⁺ T 세포 Perturb-seq atlas는 4명의 공여자로부터 얻은 약 2,200만 개의 1차 세포에서, 안정 상태(resting), 8시간, 48시간 재자극(restimulated) 조건에 걸쳐 발현 유전자와 전사인자의 합집합을 교란하였으며, 그 핵심 보고 결과는 정확히 이 프로젝트가 근거로 삼는 축이다: 활성 조절인자와 이들이 통제하는 유전자 프로그램은 활성화 맥락에 따라 상당히 달라진다([Zhu 2025, bioRxiv](https://doi.org/10.64898/2025.12.23.696273)). 이는 파이프라인에 두 가지 결과를 가져온다. 맥락-반전(context-flip) 질문 — 안정 상태와 자극 사이에서 효과 방향이 반전되는 조절인자 — 은 데이터에 내재된 것이며 atlas 자체의 핵심 결과이지 프로젝트가 임의로 부과한 것이 아니다. 그리고 배포된 `DE_stats` 객체가 이미 조건별 `log_fc`/`zscore` 레이어를 제공하므로, 프로젝트는 atlas를 즉시 사용 가능한 response-vector 행렬로 다룰 수 있으며 2,200만 세포를 재처리할 필요가 없다 — 이것이 1주일 규모의 작업 범위를 실현 가능하게 만드는 요인이다.

## 2. Perturbation-to-cell-state 투영(projection) 모델

프로젝트의 세 번째 단계 — perturbation response vector를 관찰된 조직 상태에 투영하는 것 — 는 세포 맥락 전반에 걸쳐 perturbation 효과를 예측하거나 전이(transfer)하는 것을 다루는 빠르게 발전하는 하위 분야에 속한다. 잠재(latent) 생성 모델이 관측되지 않은 맥락에서 단일세포 perturbation 반응을 예측할 수 있음을 보인 기초적 시연은 scGen이다([Lotfollahi 2019](https://doi.org/10.1038/s41592-019-0494-8)); perturbation 반응이 구조화된, 탐색 가능한 다양체(manifold) 위에 존재하여 조합과 보간(interpolation)이 의미를 갖는다는 보완적 통찰은 genetic-interaction-manifold 연구에서 나왔다([Norman 2019](https://doi.org/10.1126/science.aax4438)). 이후의 모델들은 생성 모델 접근법을 정교화하였고(적대적 생성: [Wei 2022](https://doi.org/10.1093/bioinformatics/btac357); 확산 모델: [He 2025](https://doi.org/10.1038/s41592-025-02877-y)), 특기할 만하게 *인과적(causal)* perturbation 효과를 교란된(confounded) 상관관계로부터 분리하기 시작하였으며([An 2025](https://doi.org/10.1016/j.cels.2025.101443)), 이와 병행하여 세포 간 반응의 이질성(heterogeneity)을 해석하고 해독하려는 노력도 이루어지고 있다([Song 2025](https://doi.org/10.1038/s41556-025-01626-9)). atlas 자체의 `pert2state`가 속하는 부류인 해석 가능한 state↔perturbation 대응 모델은, 이 프로젝트의 투명한 baseline과 비교되어야 할 바로 그 참조 구현이다([Xu 2024, bioRxiv](https://doi.org/10.1101/2024.03.14.585078)).

이 문헌에서 나온 두 가지 유의점이 프로젝트의 주장을 형성해야 한다. 첫째, 일반화(generalization)는 실제로 어렵고 종종 과장되어 있다: 최근의 벤치마크는 많은 perturbation-response 예측 모델이 분포 외(out-of-distribution) 상황에서 단순한 baseline을 능가하지 못함을 보여주며([Wei 2025](https://doi.org/10.1038/s41592-025-02980-0)), 커뮤니티 예측 챌린지는 정확히 이 분야가 안정적이고 정직한 성능 baseline을 갖추지 못했기 때문에 조직되었다([Zhang 2026, bioRxiv](https://doi.org/10.64898/2026.05.21.726863)). 이는 프로젝트가 투명한 dot-product 투영을 산출물로 삼고 학습된(learned) 모델은 선택적이고 명확히 표시된 비교 대상으로 다루기로 한 결정에 대한 가장 강력한 외부적 정당화이다. 둘째, 새롭게 등장하는 foundation-model 및 다층 오믹스(multilayer-omics) 접근법([Bai 2026](https://doi.org/10.1093/bioinformatics/btag307); [Chen 2026](https://doi.org/10.3389/fimmu.2026.1857447))은 유망하지만, 이 특정한 전이(dissociated in-vitro 스크린 → in-situ 종양 niche)에 대해서는 아직 입증되지 않았으며, 이를 채택하면 프로젝트가 범위 내에서 감당할 수 없는 평가 부담을 들여오게 된다. 이 계획이 이미 채택하고 있는 정직한 프레이밍 — "knockout이 CD4⁺ T 세포를 밀어붙일 전사 프로그램"이지 결코 "knockout이 niche를 리모델링할 것이다"가 아님 — 은, 이 하위 분야 자체가 자신의 전이 주장에 대해 얼마나 신중하게 유보하는지를 보면 잘 뒷받침된다.

## 3. 공간 자기상관(spatial-autocorrelation) 통계와 그 실패 양상

1단계는 공간 자기상관 — "공간상에서 비무작위적으로 조직된 발현을 보이는 유전자"에 대한 통계적 정식화 — 에 근거한다. 이 분야는 대체로 소수의 통계량(Moran's I, Geary's C, 그리고 Gaussian-process 변형)과 성숙한 구현체들로 표준화되어 있다. SpatialDE는 공간적으로 변이하는 유전자(spatially variable gene)에 대한 Gaussian-process 검정을 도입하였고([Svensson 2018](https://doi.org/10.1038/nmeth.4636)); squidpy는 순열 검정(permutation testing)과 함께 Moran's I를 패키징하였으며, 프로그램 탐지와 niche 정의 모두에 필요한 neighborhood-graph 및 niche-enrichment 기계장치를 제공한다([Palla 2022](https://doi.org/10.1038/s41592-021-01358-2)); Hotspot은 이 문제를 자기상관 기반의 *모듈(module)* 탐지로 재구성하였는데, 이는 순위가 매겨진 유전자 목록이 아니라 한 단계에서 응집된 유전자 프로그램을 산출한다는 점에서 이 프로젝트에 매력적이다([DeTomaso 2021](https://doi.org/10.1016/j.cels.2021.04.005)).

최근 방법론 문헌에서 얻는 핵심 메시지는 이러한 도구들이 **서로 일치하지 않는다**는 것이며, 이 선택은 사소한 문제가 아니다. 다수의 독립적인 2024–2025년 벤치마크는 공간적으로 변이하는 유전자 검출기(caller)들 간에, 이들이 지목하는 유전자와 위양성(false-discovery) 거동 모두에서 상당한 불일치가 있음을 보고한다([Chen 2024](https://doi.org/10.1186/s13059-023-03145-y); [Chen 2025](https://doi.org/10.1093/bioinformatics/btaf131)). 바로 이 때문에 순열(permutation) p-value를 함께 보고하는 Moran's I를 사용하고 두 번째 방법으로 교차 검증하려는 프로젝트의 계획은, 부수적인 배려가 아니라 방어 가능한 태도이다. 더 새로운 검출기들은 특정 약점을 겨냥한다: PROST와 DESpace는 검정력과 도메인 인식(domain-awareness)을 개선하고([Liang 2024](https://doi.org/10.1038/s41467-024-44835-w); [Cai 2024](https://doi.org/10.1093/bioinformatics/btae027)), 베이지안(Bayesian) 및 융합-공발현(fused-coexpression) 모델은 공간적 평활화(smoothing)와 공발현 구조를 다루며([Jiang 2022](https://doi.org/10.1002/sim.9530); [Seal 2025, bioRxiv](https://doi.org/10.1101/2025.03.29.646124)), 그리고 — 이 프로젝트의 핵심 교란 요인(confound)과 가장 관련이 깊게 — 별도의 연구 흐름은 **세포유형-특이적(cell-type-specific)** 공간 변이 유전자를 식별하여, "프로그램이 공간적으로 조직되어 있다"는 것과 "세포유형이 공간적으로 조직되어 있다"는 것을 명시적으로 구분한다([Shang 2025](https://doi.org/10.1038/s41467-025-56280-4)). 이 구분은 혼합 조직에서의 T 세포 프로그램 분석에 대한 단연 가장 중요한 방법론적 안전장치이며, 이 프로젝트가 (탈컨볼루션이 필요한 spot 단위 데이터가 아니라) 세포별 `cell_type` 라벨이 있는 단일세포 해상도 플랫폼을 사용한다는 점이 위험을 실질적으로 낮추는 이유이다. 일반적인 리뷰와 실무 가이드는 나머지 선택지들을 정리해 준다([Das 2024](https://doi.org/10.1016/j.csbj.2024.01.016); [SINFONIA, Jiang 2023](https://doi.org/10.3390/cells12040604); [Elhossiny 2026](https://doi.org/10.1007/978-1-0716-5027-1_12)).

## 4. 종양 공간 기질(substrate): NSCLC T 세포 구조와 TLS

이 다리(bridge)의 공간 쪽 절반은 조직화된 T 세포 프로그램이 실제로 존재함을 입증할 수 있는 조직을 필요로 하며, NSCLC/면역 공간 문헌은 그 기질과 내장된 양성 대조군(positive control)을 함께 제공한다. 공간적으로 해상된(spatially resolved) 연구들은 종양 내 T 세포 상태가 균일하지 않고 공간적으로 패턴화되어 있음을 확립하였다: CD8 분화와 기능장애는 공간적 위치와 niche에 의해 통제되고([Tooley 2022](https://doi.org/10.1016/j.trecan.2022.04.003); [Li 2018](https://doi.org/10.1016/j.cell.2018.11.043)), TCR-전사체 공간 지도는 뚜렷한 클론성(clonal) niche를 드러내며([Liu 2022](https://doi.org/10.1016/j.immuni.2022.09.002)), 소진(exhaustion)-연관 미세환경에 대한 단일세포/공간 지도는 면역 이웃(neighborhood)에 대한 참조 주석을 제공한다([Tietscher 2023](https://doi.org/10.1038/s41467-022-35238-w)). NSCLC에서 구체적으로는, 통합된 단일세포-및-공간 분석이 면역 미세환경과 그 골수계(myeloid) 및 림프계(lymphoid) 구획을 특성화하고([De 2024](https://doi.org/10.1038/s41467-024-48700-8); [Larroquette 2022](https://doi.org/10.1136/jitc-2021-003890); [Zhu 2025](https://doi.org/10.1016/j.ccell.2025.04.003)), 공간 시그니처는 점점 더 면역치료 결과와 연관되고 있다([Aung 2025](https://doi.org/10.1038/s41588-025-02351-7)).

"공간적으로 조직된 T 세포 프로그램"에 대한 가장 명확한 테스트 사례는 3차 림프 구조(tertiary lymphoid structure, TLS)이다 — 그 정의 자체가 공간적인, 조직화된 이소성(ectopic) 림프 집합체이다. 최근 연구는 성숙한 TLS가 종양 내(intratumoral) T 세포 및 B 세포 반응을 협조적으로 유발함을 보였고([Li 2025](https://doi.org/10.1038/s41467-025-59341-w)), 그 공간적 면역 지형을 매핑하였다([Xie 2025](https://doi.org/10.1186/s12916-025-03889-3)); 따라서 CosMx NSCLC 데이터셋에 미리 주석된 "lymphoid structure" niche는 1단계에 대한 임의의 영역이 아니라 원칙에 근거한 양성 대조군이다. 플랫폼 자체에 관해서는, 영상 기반(imaging-based) 공간 전사체학에 대한 체계적 벤치마킹이 CosMx급 어세이의 민감도와 유전자 패널 트레이드오프를 정량화하였는데([Wang 2025](https://doi.org/10.1038/s41467-025-64990-y)), 이는 프로젝트가 약 500개 유전자의 공유 공간(shared space)에서만 정렬(alignment)을 계산하고 모든 점수와 함께 그 중첩(overlap)을 보고하기로 한 결정과 직접 관련된다. 세포 이웃(cell-neighborhood) 분석 프레임워크는 투영 단계가 사용할 niche 구획화 어휘를 제공한다([Shimonov 2025, bioRxiv](https://doi.org/10.1101/2023.12.30.573739); [Lin 2023](https://doi.org/10.1016/j.cell.2022.12.028)).

## 5. 정렬(alignment), 인리치먼트(enrichment), 그리고 탈컨볼루션(deconvolution) 교란 요인

2단계의 핵심인 유사도 계산 — perturbation의 response vector가 공간 프로그램의 유전자 가중치와 얼마나 잘 일치하는지를 점수화하는 것 — 은 잘 이해된 도구를 갖춘 시그니처-인리치먼트 또는 벡터-상관 연산이다. Gene-set enrichment analysis(GSEA)는 순위가 매겨진 반응 내에서 어떤 유전자 세트가 일관되게 상향 또는 하향 조절되는지를 묻는 표준적인 방법으로 남아 있으며([Subramanian 2005](https://doi.org/10.1073/pnas.0506580102)), decoupleR와 같은 현대적 앙상블 프레임워크는 공유 유전자 공간 전반에서 인리치먼트와 벡터 기반 활성 점수화를 일관되게 실행하는 것을 용이하게 한다([Badia-i-Mompel 2022](https://doi.org/10.1093/bioadv/vbac016)). 프로젝트가 존중해야 할 방법론적 요점은 순환성(circularity)이다: 공간 프로그램을 정의하는 유전자가 perturbation 시그니처도 지배한다면, 정렬은 사소하게 부풀려진다 — 바로 이 때문에 이 계획의 순열(permutation) 및 공간-회전(spatial-rotation) 귀무(null) 모델과 (program, state)별 FDR가 이러한 상황에서의 표준적인 방어 수단이다.

가장 큰 잔여 위험은 세포유형 교란(confounding)이며, 이 지점에서 탈컨볼루션 문헌은 자원이자 경고이기도 하다. 참조 기반 탈컨볼루션 방법(cell2location, Tangram, 그리고 더 넓은 통합 생태계)은 spot당 세포유형 구성을 추정하며, "T 세포 상태"를 "T 세포 존재비(abundance)"로부터 분리하는 표준적인 방법이다([Kleshchevnikov 2022](https://doi.org/10.1038/s41587-021-01139-4); [Biancalani 2021](https://doi.org/10.1038/s41592-021-01264-7)). 그러나 독립적인 벤치마크는 이러한 방법들이 상당히 불일치하며, 어떤 단일 탈컨볼루션 또는 통합 방법도 지배적이지 않음을 보여준다([Li 2022](https://doi.org/10.1038/s41592-022-01480-9); [Yan 2023](https://doi.org/10.1093/bioinformatics/btac805); [Sang-Aram 2024](https://doi.org/10.7554/eLife.88431); [Slabowska 2024](https://doi.org/10.3389/fbinf.2024.1352594)) — 따라서 상위 후보(top hit)가 특정 탈컨볼루션 선택에 좌우되는 파이프라인은 취약하다. 이는 단일세포 해상도 경로를 택하는 것에 대한 두 번째 독립적 근거이다: 세포별 `cell_type` 라벨이 있으면 프로젝트는 CD4 memory/naive 및 Treg 세포로 직접 제한할 수 있고, 탈컨볼루션 방법 선택이라는 복불복(lottery)을 대체로 우회하여, 이 교란 요인을 모델링 문제에서 필터링 단계로 전환할 수 있다.

## 6. Perturbation과 공간(space)을 결합한 직접적 선례

두 가지 발전은 이 프로젝트의 핵심 아이디어가 특이한 것이 아니라 시의적절한 것임을 보여주며, 둘 다 가장 가까운 선행 연구(prior art)로 인용되어야 한다. 첫째, perturbation 스크리닝과 공간 readout은 이제 *실험적으로* 융합되고 있다: 동시적 CRISPR 스크리닝과 공간 전사체학은 조직 구조 내에서 perturbation이 어디에서 작용하는지를 해상하며([Binan 2025](https://doi.org/10.1016/j.cell.2025.02.012)), 대규모 병렬 in-vivo Perturb-seq은 온전한 조직에서 세포유형-특이적 perturbation 효과를 읽어낸다([Zheng 2024](https://doi.org/10.1016/j.cell.2024.04.050); 프로토콜: [Zheng 2025](https://doi.org/10.1038/s41596-024-01119-3)). 이들은 이 프로젝트의 *전산적* 다리가 저비용의 가설 생성 대리물로 기능하는 실험적 최전선을 정의하며 — 이 대조는 정직한 프레이밍을 더욱 뚜렷하게 한다: in-silico 정렬은 이후 그러한 실험이 검증할 후보의 우선순위를 매긴다. 둘째, T 세포 프로그래밍 전사인자에 대한 atlas 기반 in-silico 지명(nomination)([Chung 2025, bioRxiv](https://doi.org/10.1101/2023.01.03.522354))은 가까운 방법론적 형제이다: 이는 perturbation/발현 atlas를 사용하여 원하는 세포 상태에 대한 조절인자의 순위를 매기는데, 이는 구조적으로 이 프로젝트가 공간 프로그램으로 방향을 재설정하여 수행하는 것과 동일한 우선순위화 작업이다.

## 7. 종합 및 정직한 결론

이 파이프라인은 각 단계별로 잘 근거를 갖추고 있다: perturbation 자원은 자신의 핵심 결과가 활성화-맥락 축을 정당화하는 성숙한, 맥락-해상 atlas이며; 공간 자기상관 단계는 표준적이고 벤치마크된 통계량을 사용하고; 정렬 단계는 확립된 귀무 모델 방어 수단을 갖춘 정당한 인리치먼트/상관 계산이며; 투영 단계는 직접적인 사내(in-house) 선례와, 학습된 모델보다 투명한 baseline을 옹호하는 경계성 벤치마크 문헌을 갖고 있다. 두 가지 진정한 방법론적 위험 — **세포유형 교란**(단일세포 해상도 + 세포별 라벨 + 세포유형-특이적 SVG 인식으로 대응)과 **작은 공유 유전자 공간에서의 불안정한 정렬**(점수당 중첩을 보고하고 지지 근거가 얇은 프로그램의 가중치를 낮춤으로써 대응) — 은 모두 1차 문헌에서 반복적으로 언급되며, 모두 모델링 트릭이 아니라 프로젝트의 데이터셋 선택에 의해 완화된다. 문헌은 어떠한 산출물도 knockout이 niche를 리모델링한다는 증거로 읽는 것을 뒷받침하지 않는다 — 분리된(dissociated) 사이토카인 유도 단일배양에서 in-situ 종양으로의 전이 간극은 실재하며 아직 메워지지 않았다 — 따라서 "탐색적 투영이지, 실험적으로 검증된 것이 아니다"라는 이 계획의 고집은 상투적인 문구가 아니라, 해당 분야 자체의 유보 태도가 요구하는 입장이다. 깨끗하고 재현 가능하며 정직하게 범위가 설정된 우선순위화 도구가 달성 가능하다는 확신: 높음. 특정 상위 순위 조절인자가 실제 niche 리모델러라는 확신: 설계상 낮음 — 이것이 바로 의도된 후속 실험이며, §6에서 보았듯 이 분야는 이제 그러한 실험을 수행할 준비가 되어 있다.

---

## 참고문헌 부록

*77편의 1차 참고문헌, 방법론적 축별로 분류. DOI는 인용된 논문으로 연결된다.*


### A — Perturb-seq 및 풀링 단일세포 CRISPR 방법론

- Wang et al. 2013. *Genetic screens in human cells using the CRISPR-Cas9 system.*. Science (New York, N.Y.). [10.1126/science.1246981](https://doi.org/10.1126/science.1246981)
- Dixit et al. 2016. *Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens.*. Cell. [10.1016/j.cell.2016.11.038](https://doi.org/10.1016/j.cell.2016.11.038)
- Adamson et al. 2016. *A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response.*. Cell. [10.1016/j.cell.2016.11.048](https://doi.org/10.1016/j.cell.2016.11.048)
- Datlinger et al. 2017. *Pooled CRISPR screening with single-cell transcriptome readout.*. Nature methods. [10.1038/nmeth.4177](https://doi.org/10.1038/nmeth.4177)
- Duan et al. 2019. *Model-based understanding of single-cell CRISPR screening.*. Nature communications. [10.1038/s41467-019-10216-x](https://doi.org/10.1038/s41467-019-10216-x)
- Replogle et al. 2020. *Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing.*. Nature biotechnology. [10.1038/s41587-020-0470-y](https://doi.org/10.1038/s41587-020-0470-y)
- Replogle et al. 2022. *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq.*. Cell. [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013)
- Belk et al. 2022. *Genome-wide CRISPR screens of T cell exhaustion identify chromatin remodeling factors that limit T cell persistence.*. Cancer cell. [10.1016/j.ccell.2022.06.001](https://doi.org/10.1016/j.ccell.2022.06.001)
- Hou et al. 2022. *Single-cell CRISPR immune screens reveal immunological roles of tumor intrinsic factors.*. NAR cancer. [10.1093/narcan/zcac038](https://doi.org/10.1093/narcan/zcac038)
- Otto et al. 2023. *Structural and functional properties of mSWI/SNF chromatin remodeling complexes revealed through single-cell perturbation screens.*. Molecular cell. [10.1016/j.molcel.2023.03.013](https://doi.org/10.1016/j.molcel.2023.03.013)
- Zheng et al. 2024. *Massively parallel in vivo Perturb-seq reveals cell-type-specific transcriptional networks in cortical development.*. Cell. [10.1016/j.cell.2024.04.050](https://doi.org/10.1016/j.cell.2024.04.050)
- Pacalin et al. 2024. *Bidirectional epigenetic editing reveals hierarchies in gene regulation.*. Nature biotechnology. [10.1038/s41587-024-02213-3](https://doi.org/10.1038/s41587-024-02213-3)
- Walton et al. 2025. *CROPseq-multi: a universal solution for multiplexed perturbation in high-content pooled CRISPR screens.*. bioRxiv : the preprint server for biology. [10.1101/2024.03.17.585235](https://doi.org/10.1101/2024.03.17.585235)
- Zheng et al. 2025. *Massively parallel in vivo Perturb-seq screening.*. Nature protocols. [10.1038/s41596-024-01119-3](https://doi.org/10.1038/s41596-024-01119-3)
- Binan et al. 2025. *Simultaneous CRISPR screening and spatial transcriptomics reveal intracellular, intercellular, and functional transcriptional circuits.*. Cell. [10.1016/j.cell.2025.02.012](https://doi.org/10.1016/j.cell.2025.02.012)
- Nadig et al. 2025. *Transcriptome-wide analysis of differential expression in perturbation atlases.*. Nature genetics. [10.1038/s41588-025-02169-3](https://doi.org/10.1038/s41588-025-02169-3)
- Zhu et al. 2025. *Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits*. bioRxiv. [10.64898/2025.12.23.696273](https://doi.org/10.64898/2025.12.23.696273)

### B — CD4/T 세포 조절인자 및 종양 면역

- Shifrut et al. 2018. *Genome-wide CRISPR Screens in Primary Human T Cells Reveal Key Regulators of Immune Function.*. Cell. [10.1016/j.cell.2018.10.024](https://doi.org/10.1016/j.cell.2018.10.024)
- Li et al. 2018. *Dysfunctional CD8 T Cells Form a Proliferative, Dynamically Regulated Compartment within Human Melanoma.*. Cell. [10.1016/j.cell.2018.11.043](https://doi.org/10.1016/j.cell.2018.11.043)
- Cillo et al. 2020. *Immune Landscape of Viral- and Carcinogen-Driven Head and Neck Cancer.*. Immunity. [10.1016/j.immuni.2019.11.014](https://doi.org/10.1016/j.immuni.2019.11.014)
- Chung et al. 2025. *Atlas-Guided Discovery of Transcription Factors for T Cell Programming.*. bioRxiv : the preprint server for biology. [10.1101/2023.01.03.522354](https://doi.org/10.1101/2023.01.03.522354)
- Li et al. 2025. *Mature tertiary lymphoid structures evoke intra-tumoral T and B cell responses via progenitor exhausted CD4T cells in head and neck cancer.*. Nature communications. [10.1038/s41467-025-59341-w](https://doi.org/10.1038/s41467-025-59341-w)
- Bangs et al. 2025. *Tertiary lymphoid structures support the development of allergen-specific progenitor CD4+ T cells.*. bioRxiv : the preprint server for biology. [10.1101/2025.10.03.680350](https://doi.org/10.1101/2025.10.03.680350)

### C — Perturbation→cell-state 투영(projection) / 가상 perturbation 모델

- Lotfollahi et al. 2019. *scGen predicts single-cell perturbation responses.*. Nature methods. [10.1038/s41592-019-0494-8](https://doi.org/10.1038/s41592-019-0494-8)
- Norman et al. 2019. *Exploring genetic interaction manifolds constructed from rich single-cell phenotypes.*. Science (New York, N.Y.). [10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438)
- Wei et al. 2022. *scPreGAN, a deep generative model for predicting the response of single-cell expression to perturbation.*. Bioinformatics (Oxford, England). [10.1093/bioinformatics/btac357](https://doi.org/10.1093/bioinformatics/btac357)
- Xu et al. 2024. *Modeling interpretable correspondence between cell state and perturbation response with CellCap.*. bioRxiv : the preprint server for biology. [10.1101/2024.03.14.585078](https://doi.org/10.1101/2024.03.14.585078)
- Song et al. 2025. *Decoding heterogeneous single-cell perturbation responses.*. Nature cell biology. [10.1038/s41556-025-01626-9](https://doi.org/10.1038/s41556-025-01626-9)
- He et al. 2025. *Squidiff: predicting cellular development and responses to perturbations using a diffusion model.*. Nature methods. [10.1038/s41592-025-02877-y](https://doi.org/10.1038/s41592-025-02877-y)
- An et al. 2025. *scCausalVI disentangles single-cell perturbation responses with causality-aware generative model.*. Cell systems. [10.1016/j.cels.2025.101443](https://doi.org/10.1016/j.cels.2025.101443)
- Wei et al. 2025. *Benchmarking algorithms for generalizable single-cell perturbation response prediction.*. Nature methods. [10.1038/s41592-025-02980-0](https://doi.org/10.1038/s41592-025-02980-0)
- Li et al. 2025. *VIP-OT: Dissecting Single-Cell Biochemical State Dynamics under Perturbation via Vibrational Painting and Optimal Transport.*. bioRxiv : the preprint server for biology. [10.64898/2025.12.09.693255](https://doi.org/10.64898/2025.12.09.693255)
- Zhang et al. 2026. *A community machine learning challenge to predict the effects of gene perturbations on T cell differentiation for cancer immunotherapy.*. bioRxiv : the preprint server for biology. [10.64898/2026.05.21.726863](https://doi.org/10.64898/2026.05.21.726863)
- Chen et al. 2026. *Integratingperturbation with multilayer omics to decode regulatory networks in cancer immunity: a new frontier in precision oncology.*. Frontiers in immunology. [10.3389/fimmu.2026.1857447](https://doi.org/10.3389/fimmu.2026.1857447)
- Bai et al. 2026. *PertAdapt: unlocking single-cell foundation models for genetic perturbation prediction via condition-sensitive adaptation.*. Bioinformatics (Oxford, England). [10.1093/bioinformatics/btag307](https://doi.org/10.1093/bioinformatics/btag307)

### D — 공간 자기상관(spatial autocorrelation) 및 공간적으로 가변적인 유전자(spatially-variable-gene)/프로그램 탐지

- Svensson et al. 2018. *SpatialDE: identification of spatially variable genes.*. Nature methods. [10.1038/nmeth.4636](https://doi.org/10.1038/nmeth.4636)
- DeTomaso et al. 2021. *Hotspot identifies informative gene modules across modalities of single-cell genomics*. Cell Systems. [10.1016/j.cels.2021.04.005](https://doi.org/10.1016/j.cels.2021.04.005)
- Jiang et al. 2022. *A Bayesian modified Ising model for identifying spatially variable genes from spatial transcriptomics data.*. Statistics in medicine. [10.1002/sim.9530](https://doi.org/10.1002/sim.9530)
- Palla et al. 2022. *Squidpy: a scalable framework for spatial omics analysis*. Nature Methods. [10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2)
- Jiang et al. 2023. *SINFONIA: Scalable Identification of Spatially Variable Genes for Deciphering Spatial Domains.*. Cells. [10.3390/cells12040604](https://doi.org/10.3390/cells12040604)
- Chen et al. 2024. *Evaluating spatially variable gene detection methods for spatial transcriptomics data.*. Genome biology. [10.1186/s13059-023-03145-y](https://doi.org/10.1186/s13059-023-03145-y)
- Liang et al. 2024. *PROST: quantitative identification of spatially variable genes and domain detection in spatial transcriptomics.*. Nature communications. [10.1038/s41467-024-44835-w](https://doi.org/10.1038/s41467-024-44835-w)
- Cai et al. 2024. *DESpace: spatially variable gene detection via differential expression testing of spatial clusters.*. Bioinformatics (Oxford, England). [10.1093/bioinformatics/btae027](https://doi.org/10.1093/bioinformatics/btae027)
- Das et al. 2024. *Recent advances in spatially variable gene detection in spatial transcriptomics.*. Computational and structural biotechnology journal. [10.1016/j.csbj.2024.01.016](https://doi.org/10.1016/j.csbj.2024.01.016)
- Shang et al. 2025. *Statistical identification of cell type-specific spatially variable genes in spatial transcriptomics.*. Nature communications. [10.1038/s41467-025-56280-4](https://doi.org/10.1038/s41467-025-56280-4)
- Chen et al. 2025. *Benchmarking algorithms for spatially variable gene identification in spatial transcriptomics.*. Bioinformatics (Oxford, England). [10.1093/bioinformatics/btaf131](https://doi.org/10.1093/bioinformatics/btaf131)
- Seal et al. 2025. *SpaceBF: Spatial coexpression analysis using Bayesian Fused approaches in spatial omics datasets.*. bioRxiv : the preprint server for biology. [10.1101/2025.03.29.646124](https://doi.org/10.1101/2025.03.29.646124)
- Elhossiny et al. 2026. *From Raw Data to Biological Insights: A Practical Guide for Spatial Transcriptomics Analysis in R and Python.*. Methods in molecular biology (Clifton, N.J.). [10.1007/978-1-0716-5027-1_12](https://doi.org/10.1007/978-1-0716-5027-1_12)

### E — 면역 구조에 대한 종양 공간 전사체학(NSCLC, TLS)

- Tooley et al. 2022. *Spatial determinants of CD8T cell differentiation in cancer.*. Trends in cancer. [10.1016/j.trecan.2022.04.003](https://doi.org/10.1016/j.trecan.2022.04.003)
- Larroquette et al. 2022. *Spatial transcriptomics of macrophage infiltration in non-small cell lung cancer reveals determinants of sensitivity and resistance to anti-PD1/PD-L1 antibodies.*. Journal for immunotherapy of cancer. [10.1136/jitc-2021-003890](https://doi.org/10.1136/jitc-2021-003890)
- Liu et al. 2022. *Spatial maps of T cell receptors and transcriptomes reveal distinct immune niches and interactions in the adaptive immune response.*. Immunity. [10.1016/j.immuni.2022.09.002](https://doi.org/10.1016/j.immuni.2022.09.002)
- Tietscher et al. 2023. *A comprehensive single-cell map of T cell exhaustion-associated immune environments in human breast cancer.*. Nature communications. [10.1038/s41467-022-35238-w](https://doi.org/10.1038/s41467-022-35238-w)
- Lin et al. 2023. *Multiplexed 3D atlas of state transitions and immune interaction in colorectal cancer.*. Cell. [10.1016/j.cell.2022.12.028](https://doi.org/10.1016/j.cell.2022.12.028)
- Hashimoto et al. 2024. *Spatial and single-cell colocalisation analysis reveals MDK-mediated immunosuppressive environment with regulatory T cells in colorectal carcinogenesis.*. EBioMedicine. [10.1016/j.ebiom.2024.105102](https://doi.org/10.1016/j.ebiom.2024.105102)
- De et al. 2024. *Single-cell and spatial transcriptomics analysis of non-small cell lung cancer.*. Nature communications. [10.1038/s41467-024-48700-8](https://doi.org/10.1038/s41467-024-48700-8)
- Peyraud et al. 2025. *Spatially resolved transcriptomics reveal the determinants of primary resistance to immunotherapy in NSCLC with mature tertiary lymphoid structures.*. Cell reports. Medicine. [10.1016/j.xcrm.2025.101934](https://doi.org/10.1016/j.xcrm.2025.101934)
- Shimonov et al. 2025. *SORBET: Automated cell-neighborhood analysis of spatial transcriptomics or proteomics for interpretable sample classification via GNN.*. bioRxiv : the preprint server for biology. [10.1101/2023.12.30.573739](https://doi.org/10.1101/2023.12.30.573739)
- Xie et al. 2025. *Immune microenvironment spatial landscapes of tertiary lymphoid structures in gastric cancer.*. BMC medicine. [10.1186/s12916-025-03889-3](https://doi.org/10.1186/s12916-025-03889-3)
- Pentimalli et al. 2025. *Combining spatial transcriptomics and ECM imaging in 3D for mapping cellular interactions in the tumor microenvironment.*. Cell systems. [10.1016/j.cels.2025.101261](https://doi.org/10.1016/j.cels.2025.101261)
- Zhu et al. 2025. *Spatial and multiomics analysis of human and mouse lung adenocarcinoma precursors reveals TIM-3 as a putative target for precancer interception.*. Cancer cell. [10.1016/j.ccell.2025.04.003](https://doi.org/10.1016/j.ccell.2025.04.003)
- Aung et al. 2025. *Spatial signatures for predicting immunotherapy outcomes using multi-omics in non-small cell lung cancer.*. Nature genetics. [10.1038/s41588-025-02351-7](https://doi.org/10.1038/s41588-025-02351-7)
- Cakmak et al. 2025. *Spatial immune profiling defines a subset of human gliomas with functional tertiary lymphoid structures.*. Immunity. [10.1016/j.immuni.2025.09.018](https://doi.org/10.1016/j.immuni.2025.09.018)
- Wang et al. 2025. *Systematic benchmarking of imaging spatial transcriptomics platforms in FFPE tissues.*. Nature communications. [10.1038/s41467-025-64990-y](https://doi.org/10.1038/s41467-025-64990-y)
- Zhai et al. 2025. *Spatial proteomic profiling reveals conserved prognostic immune microenvironment features across molecular subtypes in small cell lung cancer.*. Pharmacological research. [10.1016/j.phrs.2025.108048](https://doi.org/10.1016/j.phrs.2025.108048)
- Wang et al. 2025. *Single-cell and spatial transcriptomics implicate a prognostic function of tertiary lymphoid structures in gastric cancer.*. Nature communications. [10.1038/s41467-025-65421-8](https://doi.org/10.1038/s41467-025-65421-8)
- Wang et al. 2026. *Integrative Multi-Omics Analysis Identifies an-Associated Spatial Mesenchymal-Myeloid Program in Glioblastoma.*. Genes. [10.3390/genes17060610](https://doi.org/10.3390/genes17060610)

### F — 시그니처 정렬, enrichment 및 deconvolution 방법

- Subramanian et al. 2005. *Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles*. PNAS. [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
- Biancalani et al. 2021. *Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram*. Nature Methods. [10.1038/s41592-021-01264-7](https://doi.org/10.1038/s41592-021-01264-7)
- Li et al. 2022. *Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution.*. Nature methods. [10.1038/s41592-022-01480-9](https://doi.org/10.1038/s41592-022-01480-9)
- Kleshchevnikov et al. 2022. *Cell2location maps fine-grained cell types in spatial transcriptomics*. Nature Biotechnology. [10.1038/s41587-021-01139-4](https://doi.org/10.1038/s41587-021-01139-4)
- Badia-i-Mompel et al. 2022. *decoupleR: ensemble of computational methods to infer biological activities from omics data*. Bioinformatics Advances. [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)
- Yan et al. 2023. *Benchmarking and integration of methods for deconvoluting spatial transcriptomic data.*. Bioinformatics (Oxford, England). [10.1093/bioinformatics/btac805](https://doi.org/10.1093/bioinformatics/btac805)
- Slabowska et al. 2024. *A systematic evaluation of state-of-the-art deconvolution methods in spatial transcriptomics: insights from cardiovascular disease and chronic kidney disease.*. Frontiers in bioinformatics. [10.3389/fbinf.2024.1352594](https://doi.org/10.3389/fbinf.2024.1352594)
- Sang-Aram et al. 2024. *Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics.*. eLife. [10.7554/eLife.88431](https://doi.org/10.7554/eLife.88431)
- Xi et al. 2025. *SegDecon bridges histology and transcriptomics through AI-based nuclei segmentation and image-informed spatial deconvolution.*. Computational and structural biotechnology journal. [10.1016/j.csbj.2025.10.041](https://doi.org/10.1016/j.csbj.2025.10.041)
- Zhang et al. 2025. *Single-cell/spatial integration reveals an MES2-like glioblastoma program orchestrated by immune communication and regulatory networks.*. Frontiers in immunology. [10.3389/fimmu.2025.1699134](https://doi.org/10.3389/fimmu.2025.1699134)
- Wu et al. 2026. *Integrative single-cell and spatial transcriptomics with machine learning identify a Luminal-inflam malignant program and reveal an RPN1-PERK UPR vulnerability in triple-negative breast cancer.*. Scientific reports. [10.1038/s41598-026-56434-4](https://doi.org/10.1038/s41598-026-56434-4)
