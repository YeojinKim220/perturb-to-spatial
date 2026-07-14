import React from "react";
import { interpolate, useCurrentFrame } from "remotion";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

export const Title: React.FC = () => {
  const frame = useCurrentFrame();
  const a = useReveal(6);
  const b = useReveal(20);
  const c = useReveal(40);
  const line = interpolate(frame, [10, 45], [0, 1], {
    extrapolateLeft: "clamp",
    extrapolateRight: "clamp",
  });
  return (
    <Bg>
      <div
        style={{
          position: "absolute",
          inset: 0,
          display: "flex",
          flexDirection: "column",
          justifyContent: "center",
          padding: "0 140px",
        }}
      >
        <div style={a}>
          <Kicker>Hackathon Report · Researcher Track</Kicker>
        </div>
        <div
          style={{
            ...b,
            marginTop: 34,
            fontSize: 66,
            fontWeight: 800,
            lineHeight: 1.14,
            letterSpacing: -1,
            maxWidth: 1500,
          }}
        >
          Do naive-annotated CD4⁺ T cells occupy distinct{" "}
          <span style={{ color: theme.accent }}>spatial contexts</span> as
          distinct <span style={{ color: theme.accent2 }}>functional states</span>{" "}
          in NSCLC?
        </div>
        <div
          style={{
            marginTop: 40,
            height: 3,
            width: `${line * 620}px`,
            background: `linear-gradient(90deg, ${theme.accent}, ${theme.accent2})`,
            borderRadius: 2,
          }}
        />
        <div
          style={{
            ...c,
            marginTop: 40,
            fontSize: 32,
            color: theme.sub,
            maxWidth: 1350,
            lineHeight: 1.5,
          }}
        >
          Building a transparent bridge from a genome-scale CD4⁺ Perturb-seq atlas
          to tumor spatial transcriptomics — and asking whether it can nominate
          regulators for each state.
        </div>
      </div>
      <div
        style={{
          position: "absolute",
          bottom: 70,
          left: 140,
          fontFamily: theme.fontMono,
          fontSize: 24,
          color: theme.faint,
        }}
      >
        human · GRCh38 · CosMx NSCLC × Marson–Pritchard CD4⁺ Perturb-seq
      </div>
    </Bg>
  );
};
