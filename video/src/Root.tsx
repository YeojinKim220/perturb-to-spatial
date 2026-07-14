import React from "react";
import { Composition, Series, AbsoluteFill } from "remotion";
import { theme, fps } from "./theme";
import { Title } from "./scenes/Title";
import { Datasets } from "./scenes/Datasets";
import { FigureScene } from "./scenes/FigureScene";
import { WideFigureScene } from "./scenes/WideFigureScene";
import { Funnel } from "./scenes/Funnel";
import { Conclusion } from "./scenes/Conclusion";
import { Closing } from "./scenes/Closing";

// durations in frames @ 30fps
const D = {
  title: 165,
  datasets: 255,
  fig0a: 285,
  fig0b: 270,
  fig1a: 285,
  fig1: 345,
  fig2a: 315,
  fig2b: 300,
  fig3: 315,
  funnel: 285,
  conclusion: 330,
  closing: 255,
};

const TOTAL = Object.values(D).reduce((a, b) => a + b, 0); // 3405f ≈ 113.5s

const Report: React.FC = () => (
  <AbsoluteFill style={{ backgroundColor: theme.bg }}>
    <Series>
      <Series.Sequence durationInFrames={D.title}>
        <Title />
      </Series.Sequence>
      <Series.Sequence durationInFrames={D.datasets}>
        <Datasets />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig0a}>
        <WideFigureScene
          index="Fig 0a"
          kicker="The spatial substrate"
          title="Three LUAD-5 tumor sections, 22 cell types"
          img="fig0a_celltype_1x3.png"
          accent={theme.accent}
          chips={[
            { text: "765,771 cells", color: theme.accent },
            { text: "960-gene panel", color: theme.accent2 },
            { text: "CosMx NSCLC", color: theme.gold },
          ]}
          caption="Every cell colored by its annotated type — tumor nests, stromal bands, and immune infiltrate. These are the tissue maps every downstream spatial analysis is computed on."
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig0b}>
        <WideFigureScene
          index="Fig 0b"
          kicker="Focusing on the CD4⁺ T cells"
          title="Naive, memory & Treg cells thread through the tissue"
          img="fig0b_cd4_only_1x3.png"
          accent={theme.accent2}
          chips={[
            { text: "naive", color: theme.accent2 },
            { text: "memory", color: theme.gold },
            { text: "Treg", color: theme.good },
          ]}
          caption="The same sections with only CD4⁺ T-cell subtypes colored; all else grey. The <b>23,889 naive-annotated CD4⁺ cells</b> (blue) are the substrate for the state × context analysis."
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig1a}>
        <WideFigureScene
          index="Fig 1a"
          kicker="Same annotation, different context"
          title="Where the intrinsic states sit relative to tumor"
          img="fig1a_state_contact_map_1x3.png"
          accent={theme.accent}
          chips={[
            { text: "Th1-high", color: theme.warn },
            { text: "cytotoxic-high", color: theme.accent2 },
            { text: "tumor (grey)", color: theme.faint },
          ]}
          caption="Naive CD4⁺ cells colored by dominant state over a grey tumor background. A tissue-level view only — the contact-enrichment effect itself is the quantitative log-OR test in Fig 1. <i>Exploratory — not experimentally validated.</i>"
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig1}>
        <FigureScene
          index="Fig 1"
          kicker="Quantifying state × context"
          title="Weak, mostly non-reproducing spatial structure"
          img="fig1_state_context.png"
          accent={theme.accent}
          bullets={[
            {
              text: "Only the <b>tumor-distance axis</b> is informative — Th1-high enriched at contact (+0.28), cytotoxic-high depleted (−0.86), both sign-consistent 3/3.",
              color: theme.accent,
              delay: 34,
            },
            {
              text: "After graph correction, spatial autocorrelation is <b>not significant: 0 / 15</b> Moran's I tests (BH).",
              color: theme.warn,
              delay: 58,
            },
            {
              text: "Leave-one-section-out reproducibility is <b>low (ρ ≈ 0.31)</b> — the earlier 0.89–0.93 was inflated by joint niche definition.",
              color: theme.gold,
              delay: 82,
            },
          ]}
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig2a}>
        <FigureScene
          index="Fig 2a"
          kicker="Perturbation depends on activation state"
          title="Knockdown effects flip with Rest / Stim"
          img="fig2a_context.png"
          accent={theme.accent2}
          bullets={[
            {
              text: "Screening every (target, program) pair, <b>5,736 of 55,510</b> reverse direction between Rest and Stim 48 h.",
              color: theme.accent2,
              delay: 34,
            },
            {
              text: "But that count is a <b>ceiling</b>: 2,236 flips are incomplete (no significant downstream gene). Reliable count is far smaller.",
              color: theme.warn,
              delay: 58,
            },
            {
              text: "TCR-signaling knockdowns (ZAP70, CD3E, VAV1) act mainly when stimulated — a sane sanity signal.",
              color: theme.good,
              delay: 82,
            },
          ]}
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig2b}>
        <FigureScene
          index="Fig 2b"
          kicker="Per-gene rest-vs-stim response"
          title="Some regulators are stable, others context-dependent"
          img="fig2b_flips.png"
          accent={theme.accent2}
          bullets={[
            {
              text: "Each point is one downstream gene; <b>red points reverse direction</b> between rest and stim.",
              color: theme.warn,
              delay: 34,
            },
            {
              text: "<b>TP53 (r = 0.57)</b> is mostly state-stable — the diagonal holds across conditions.",
              color: theme.good,
              delay: 58,
            },
            {
              text: "<b>GATA3 (r = 0.36)</b> and <b>ZAP70 (r = 0.27)</b> grow progressively more context-dependent.",
              color: theme.accent2,
              delay: 82,
            },
          ]}
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.fig3}>
        <FigureScene
          index="Fig 3"
          kicker="Linking spatial contrast to perturbation"
          title="Stress-corrected — but not biologically specific"
          img="fig3_integration.png"
          accent={theme.gold}
          bullets={[
            {
              text: "Corrected contact-vs-distal signatures are <b>no longer stress-dominated</b> (CD74, PECAM1, ICAM1, GZMA survive).",
              color: theme.good,
              delay: 34,
            },
            {
              text: "Yet removing stress barely moves the ranking (ρ = 0.97–0.99) — aligned top targets are generic replication genes.",
              color: theme.warn,
              delay: 60,
            },
          ]}
        />
      </Series.Sequence>

      <Series.Sequence durationInFrames={D.funnel}>
        <Funnel />
      </Series.Sequence>
      <Series.Sequence durationInFrames={D.conclusion}>
        <Conclusion />
      </Series.Sequence>
      <Series.Sequence durationInFrames={D.closing}>
        <Closing />
      </Series.Sequence>
    </Series>
  </AbsoluteFill>
);

export const RemotionRoot: React.FC = () => (
  <Composition
    id="Report"
    component={Report}
    durationInFrames={TOTAL}
    fps={fps}
    width={1920}
    height={1080}
  />
);
