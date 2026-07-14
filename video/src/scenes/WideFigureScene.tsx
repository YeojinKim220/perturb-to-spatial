import React from "react";
import { Img, staticFile, interpolate, useCurrentFrame } from "remotion";
import { Bg, Kicker, useReveal } from "../ui";
import { theme } from "../theme";

// For panoramic 1x3 tissue maps: full-width image on top, caption chips below.
export const WideFigureScene: React.FC<{
  index: string;
  kicker: string;
  title: string;
  img: string;
  caption: string;
  chips?: { text: string; color: string }[];
  accent?: string;
}> = ({ index, kicker, title, img, caption, chips = [], accent = theme.accent }) => {
  const head = useReveal(4);
  const cap = useReveal(30);
  const frame = useCurrentFrame();
  const imgIn = interpolate(frame, [8, 32], [0, 1], {
    extrapolateLeft: "clamp",
    extrapolateRight: "clamp",
  });
  return (
    <Bg>
      <div style={{ padding: "58px 96px", height: "100%", boxSizing: "border-box" }}>
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
            marginTop: 40,
            opacity: imgIn,
            transform: `translateY(${interpolate(imgIn, [0, 1], [18, 0])}px)`,
            background: "#ffffff",
            borderRadius: 16,
            padding: 20,
            border: `1px solid ${theme.line}`,
            boxShadow: "0 26px 64px rgba(0,0,0,0.5)",
          }}
        >
          <Img
            src={staticFile(img)}
            style={{ width: "100%", display: "block", borderRadius: 6 }}
          />
        </div>

        <div style={{ ...cap, marginTop: 40 }}>
          <div style={{ display: "flex", gap: 14, flexWrap: "wrap", marginBottom: 20 }}>
            {chips.map((c, i) => (
              <span
                key={i}
                style={{
                  display: "inline-flex",
                  alignItems: "center",
                  gap: 10,
                  padding: "8px 18px",
                  borderRadius: 999,
                  border: `1px solid ${c.color}55`,
                  background: `${c.color}14`,
                  color: c.color,
                  fontSize: 26,
                  fontWeight: 700,
                }}
              >
                <span
                  style={{
                    width: 12,
                    height: 12,
                    borderRadius: 4,
                    background: c.color,
                    boxShadow: `0 0 14px ${c.color}88`,
                  }}
                />
                {c.text}
              </span>
            ))}
          </div>
          <div
            style={{ fontSize: 31, lineHeight: 1.45, color: theme.sub, maxWidth: 1500 }}
            dangerouslySetInnerHTML={{ __html: caption }}
          />
        </div>
      </div>
    </Bg>
  );
};
