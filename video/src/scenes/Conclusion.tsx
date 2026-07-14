import React from "react";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

const Col: React.FC<{
  title: string;
  color: string;
  items: string[];
  baseDelay: number;
}> = ({ title, color, items, baseDelay }) => {
  const h = useReveal(baseDelay);
  return (
    <div
      style={{
        flex: 1,
        background: theme.bgPanel,
        border: `1px solid ${theme.line}`,
        borderLeft: `4px solid ${color}`,
        borderRadius: 16,
        padding: "32px 36px",
      }}
    >
      <div style={{ ...h, fontSize: 30, fontWeight: 800, color }}>{title}</div>
      <div style={{ marginTop: 22, display: "flex", flexDirection: "column", gap: 18 }}>
        {items.map((it, i) => (
          <Item key={i} text={it} delay={baseDelay + 8 + i * 7} color={color} />
        ))}
      </div>
    </div>
  );
};

const Item: React.FC<{ text: string; delay: number; color: string }> = ({
  text,
  delay,
  color,
}) => {
  const r = useReveal(delay);
  return (
    <div style={{ ...r, display: "flex", gap: 14 }}>
      <div style={{ color, fontSize: 26, lineHeight: 1.35 }}>›</div>
      <div style={{ fontSize: 26, lineHeight: 1.4, color: theme.ink }}>{text}</div>
    </div>
  );
};

export const Conclusion: React.FC = () => {
  const head = useReveal(4);
  const foot = useReveal(120);
  return (
    <Bg>
      <div style={{ padding: "68px 110px" }}>
        <div style={head}>
          <Kicker color={theme.gold}>Conclusion</Kicker>
          <div style={{ fontSize: 48, fontWeight: 800, marginTop: 14, maxWidth: 1500 }}>
            The bridge works — but it selected no defensible target
          </div>
        </div>

        <div style={{ display: "flex", gap: 30, marginTop: 44 }}>
          <Col
            title="What is established"
            color={theme.good}
            baseDelay={20}
            items={[
              "A transparent, reusable pipeline: independent state × context, activation-aware screen, audited FDR + reliability + viability.",
              "Cytotoxic enhancement is the only endpoint carrying a transferable directional signal.",
              "Th1-high enriched at tumor contact; cytotoxic-high depleted (sign-consistent 3/3).",
            ]}
          />
          <Col
            title="What it is not"
            color={theme.warn}
            baseDelay={50}
            items={[
              "No endpoint yields a clean, biologically credible target shortlist.",
              "Corrected Moran's I: 0 / 15 significant — no cohesive spatial territories.",
              "Superseded claims (PTPN2, MBTPS1, CDK9) came from a pooled-null bug and are retracted.",
            ]}
          />
        </div>

        <div
          style={{
            ...foot,
            marginTop: 40,
            fontSize: 32,
            textAlign: "center",
            color: theme.ink,
            lineHeight: 1.5,
          }}
        >
          A genuine, calibrated{" "}
          <span style={{ color: theme.gold, fontWeight: 800 }}>negative result</span>:
          FDR alone cannot pick a target here. Every map is exploratory —{" "}
          <span style={{ color: theme.sub }}>not experimentally validated.</span>
        </div>
      </div>
    </Bg>
  );
};
