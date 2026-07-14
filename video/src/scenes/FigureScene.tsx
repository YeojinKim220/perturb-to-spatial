import React from "react";
import { Img, staticFile, interpolate, useCurrentFrame } from "remotion";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

export type Bullet = { text: string; color?: string; delay: number };

export const FigureScene: React.FC<{
  index: string;
  kicker: string;
  title: string;
  img: string;
  bullets: Bullet[];
  accent?: string;
}> = ({ index, kicker, title, img, bullets, accent = theme.accent }) => {
  const head = useReveal(4);
  const frame = useCurrentFrame();
  const imgIn = interpolate(frame, [8, 30], [0, 1], {
    extrapolateLeft: "clamp",
    extrapolateRight: "clamp",
  });
  return (
    <Bg>
      <div style={{ padding: "62px 96px", height: "100%", boxSizing: "border-box" }}>
        <div style={{ ...head, display: "flex", alignItems: "baseline", gap: 22 }}>
          <div
            style={{
              fontFamily: theme.fontMono,
              fontSize: 30,
              fontWeight: 800,
              color: accent,
              border: `2px solid ${accent}66`,
              borderRadius: 12,
              padding: "4px 16px",
            }}
          >
            {index}
          </div>
          <div>
            <Kicker color={accent}>{kicker}</Kicker>
            <div style={{ fontSize: 42, fontWeight: 800, marginTop: 8 }}>{title}</div>
          </div>
        </div>

        <div
          style={{
            display: "flex",
            gap: 46,
            marginTop: 34,
            alignItems: "center",
          }}
        >
          <div
            style={{
              flex: "1.35",
              opacity: imgIn,
              transform: `translateY(${interpolate(imgIn, [0, 1], [16, 0])}px)`,
              background: "#ffffff",
              borderRadius: 16,
              padding: 18,
              border: `1px solid ${theme.line}`,
              boxShadow: "0 24px 60px rgba(0,0,0,0.45)",
            }}
          >
            <Img
              src={staticFile(img)}
              style={{ width: "100%", display: "block", borderRadius: 6 }}
            />
          </div>

          <div style={{ flex: "1", display: "flex", flexDirection: "column", gap: 24 }}>
            {bullets.map((b, i) => (
              <BulletRow key={i} b={b} />
            ))}
          </div>
        </div>
      </div>
    </Bg>
  );
};

const BulletRow: React.FC<{ b: Bullet }> = ({ b }) => {
  const r = useReveal(b.delay);
  const color = b.color ?? theme.accent;
  return (
    <div style={{ ...r, display: "flex", gap: 18 }}>
      <div
        style={{
          marginTop: 12,
          minWidth: 12,
          width: 12,
          height: 12,
          borderRadius: 4,
          background: color,
          boxShadow: `0 0 16px ${color}88`,
        }}
      />
      <div
        style={{
          fontSize: 30,
          lineHeight: 1.42,
          color: theme.ink,
        }}
        dangerouslySetInnerHTML={{ __html: b.text }}
      />
    </div>
  );
};
