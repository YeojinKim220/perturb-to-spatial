import React from "react";
import { interpolate, useCurrentFrame } from "remotion";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

const items = [
  "Autonomously drove a remote HPC / Slurm cluster — built envs, submitted jobs, collected outputs.",
  "Audited its own statistics and retracted a headline claim (the impossible “FDR 0.024”).",
  "Fixed the spatial-graph & reproducibility design (Moran cutoffs; leave-one-section-out).",
  "Reported the calibrated negative result — and refused to build a misleading Figure 4.",
];

export const Closing: React.FC = () => {
  const frame = useCurrentFrame();
  const head = useReveal(4);
  const glow = interpolate(frame, [0, 60], [0, 1], { extrapolateRight: "clamp" });
  return (
    <Bg>
      <div style={{ padding: "80px 130px" }}>
        <div style={head}>
          <Kicker color={theme.accent2}>How Claude was used</Kicker>
          <div style={{ fontSize: 46, fontWeight: 800, marginTop: 14 }}>
            Human-in-the-loop, honest to the end
          </div>
        </div>
        <div style={{ marginTop: 46, display: "flex", flexDirection: "column", gap: 26 }}>
          {items.map((it, i) => (
            <Row key={i} text={it} delay={18 + i * 12} />
          ))}
        </div>
        <div
          style={{
            marginTop: 64,
            paddingTop: 34,
            borderTop: `1px solid ${theme.line}`,
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
            opacity: glow,
          }}
        >
          <div style={{ fontSize: 30, color: theme.sub }}>
            The most useful output this round is a{" "}
            <span style={{ color: theme.gold, fontWeight: 700 }}>
              calibrated negative result.
            </span>
          </div>
          <div
            style={{
              fontFamily: theme.fontMono,
              fontSize: 24,
              color: theme.faint,
            }}
          >
            2026-07-13 · human · GRCh38
          </div>
        </div>
      </div>
    </Bg>
  );
};

const Row: React.FC<{ text: string; delay: number }> = ({ text, delay }) => {
  const r = useReveal(delay);
  return (
    <div style={{ ...r, display: "flex", gap: 18, alignItems: "flex-start" }}>
      <div
        style={{
          minWidth: 14,
          width: 14,
          height: 14,
          marginTop: 10,
          borderRadius: 4,
          background: theme.accent2,
          boxShadow: `0 0 14px ${theme.accent2}88`,
        }}
      />
      <div style={{ fontSize: 31, lineHeight: 1.4, color: theme.ink, maxWidth: 1400 }}>
        {text}
      </div>
    </div>
  );
};
