import React from "react";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

const Card: React.FC<{
  title: string;
  color: string;
  rows: [string, string][];
  style?: React.CSSProperties;
}> = ({ title, color, rows, style }) => (
  <div
    style={{
      flex: 1,
      background: theme.bgPanel,
      border: `1px solid ${theme.line}`,
      borderTop: `4px solid ${color}`,
      borderRadius: 18,
      padding: "36px 40px",
      ...style,
    }}
  >
    <div style={{ fontSize: 34, fontWeight: 800, color }}>{title}</div>
    <div style={{ marginTop: 24 }}>
      {rows.map(([k, v]) => (
        <div key={k} style={{ marginBottom: 20 }}>
          <div
            style={{
              fontFamily: theme.fontMono,
              fontSize: 22,
              color: theme.faint,
              textTransform: "uppercase",
              letterSpacing: 1,
            }}
          >
            {k}
          </div>
          <div style={{ fontSize: 28, color: theme.ink, marginTop: 4, lineHeight: 1.35 }}>
            {v}
          </div>
        </div>
      ))}
    </div>
  </div>
);

export const Datasets: React.FC = () => {
  const head = useReveal(4);
  const l = useReveal(16);
  const r = useReveal(26);
  const foot = useReveal(44);
  return (
    <Bg>
      <div style={{ padding: "80px 120px" }}>
        <div style={head}>
          <Kicker>Two atlases · one shared gene space</Kicker>
          <div style={{ fontSize: 46, fontWeight: 800, marginTop: 16 }}>
            One annotated identity, joined across platforms
          </div>
        </div>
        <div style={{ display: "flex", gap: 34, marginTop: 46 }}>
          <div style={{ ...l, display: "flex", flex: 1 }}>
            <Card
              title="Perturbation atlas"
              color={theme.accent2}
              rows={[
                ["Source", "Genome-scale CD4⁺ Perturb-seq (CRISPRi, 4 donors)"],
                ["Measures", "logFC / z-score of every gene after each knockdown"],
                ["Key axis", "Rest · Stim 8 h · Stim 48 h"],
                ["Used here", "DE_stats.h5ad — 10,282 genes"],
              ]}
            />
          </div>
          <div style={{ ...r, display: "flex", flex: 1 }}>
            <Card
              title="Spatial atlas"
              color={theme.accent}
              rows={[
                ["Source", "CosMx NSCLC (He 2022)"],
                ["Measures", "Single-cell spatial location + expression"],
                ["Scale", "765,771 cells · 960-gene panel · 22 cell types"],
                ["Used here", "23,889 naive CD4⁺ cells · LUAD-5 R1–R3"],
              ]}
            />
          </div>
        </div>
        <div
          style={{
            ...foot,
            marginTop: 40,
            fontSize: 28,
            color: theme.sub,
            textAlign: "center",
            lineHeight: 1.5,
          }}
        >
          Cross-platform overlap is small and always reported:{" "}
          <span style={{ color: theme.gold, fontWeight: 700 }}>501 genes</span>{" "}
          (gene-level) ·{" "}
          <span style={{ color: theme.gold, fontWeight: 700 }}>912 genes</span>{" "}
          (state-matched). A one-week join with a strict inference ceiling.
        </div>
      </div>
    </Bg>
  );
};
