import React from "react";
import { interpolate, spring, useCurrentFrame, useVideoConfig } from "remotion";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

const steps = [
  { label: "All targets", value: "11,287", note: "genome-scale screen" },
  { label: "FDR < 0.1 (per-target)", value: "127", note: "cytotoxic endpoint" },
  { label: "Donor-reliable", value: "12", note: "reproduces across donors" },
  { label: "No viability flag", value: "3", note: "COX6B1 · ITPA · WARS2" },
];

export const Funnel: React.FC = () => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();
  const head = useReveal(4);
  const foot = useReveal(150);
  return (
    <Bg>
      <div style={{ padding: "70px 120px" }}>
        <div style={head}>
          <Kicker color={theme.warn}>Figure 3 · The evidence funnel</Kicker>
          <div style={{ fontSize: 46, fontWeight: 800, marginTop: 14 }}>
            FDR + reliability + viability, applied together
          </div>
        </div>

        <div
          style={{
            marginTop: 56,
            display: "flex",
            flexDirection: "column",
            gap: 22,
            alignItems: "center",
          }}
        >
          {steps.map((s, i) => {
            const delay = 20 + i * 26;
            const sp = spring({
              frame: frame - delay,
              fps,
              config: { damping: 200 },
              durationInFrames: 22,
            });
            const width = interpolate(i, [0, steps.length - 1], [1180, 360]);
            const isLast = i === steps.length - 1;
            const color = isLast ? theme.gold : theme.accent;
            return (
              <div
                key={i}
                style={{
                  opacity: sp,
                  transform: `scale(${interpolate(sp, [0, 1], [0.9, 1])})`,
                  width,
                  background: `linear-gradient(90deg, ${color}22, ${color}0c)`,
                  border: `1px solid ${color}66`,
                  borderRadius: 14,
                  padding: "18px 34px",
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "space-between",
                }}
              >
                <div style={{ fontSize: 27, color: theme.ink, fontWeight: 600 }}>
                  {s.label}
                </div>
                <div style={{ display: "flex", alignItems: "baseline", gap: 16 }}>
                  <div
                    style={{
                      fontFamily: theme.fontMono,
                      fontSize: 38,
                      fontWeight: 800,
                      color,
                    }}
                  >
                    {s.value}
                  </div>
                </div>
              </div>
            );
          })}
        </div>

        <div
          style={{
            ...foot,
            marginTop: 46,
            fontSize: 30,
            textAlign: "center",
            color: theme.sub,
            lineHeight: 1.5,
          }}
        >
          The 3 survivors —{" "}
          <span style={{ color: theme.gold, fontFamily: theme.fontMono }}>
            COX6B1, ITPA, WARS2
          </span>{" "}
          — are mitochondrial / tRNA / nucleotide genes with{" "}
          <span style={{ color: theme.warn, fontWeight: 700 }}>
            no cytotoxic rationale
          </span>
          . The funnel ends with no clean shortlist.
        </div>
      </div>
    </Bg>
  );
};
