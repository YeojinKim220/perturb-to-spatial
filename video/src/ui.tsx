import React from "react";
import {
  interpolate,
  spring,
  useCurrentFrame,
  useVideoConfig,
  AbsoluteFill,
} from "remotion";
import { theme } from "./theme";

export const Bg: React.FC<{ children: React.ReactNode }> = ({ children }) => (
  <AbsoluteFill
    style={{
      backgroundColor: theme.bg,
      backgroundImage:
        "radial-gradient(1100px 700px at 82% -8%, rgba(77,208,225,0.10), transparent 60%), radial-gradient(900px 600px at 8% 108%, rgba(124,156,255,0.10), transparent 60%)",
      fontFamily: theme.fontSans,
      color: theme.ink,
    }}
  >
    {children}
  </AbsoluteFill>
);

// fade + rise in, with optional out
export const useReveal = (delay = 0, dur = 20) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();
  const s = spring({
    frame: frame - delay,
    fps,
    config: { damping: 200, mass: 0.7 },
    durationInFrames: dur,
  });
  return {
    opacity: s,
    transform: `translateY(${interpolate(s, [0, 1], [18, 0])}px)`,
  };
};

export const useFadeOut = (start: number, dur = 15) => {
  const frame = useCurrentFrame();
  return interpolate(frame, [start, start + dur], [1, 0], {
    extrapolateLeft: "clamp",
    extrapolateRight: "clamp",
  });
};

export const Kicker: React.FC<{ children: React.ReactNode; color?: string }> = ({
  children,
  color = theme.accent,
}) => (
  <div
    style={{
      fontFamily: theme.fontMono,
      fontSize: 26,
      letterSpacing: 4,
      textTransform: "uppercase",
      color,
      fontWeight: 600,
    }}
  >
    {children}
  </div>
);

export const Chip: React.FC<{
  children: React.ReactNode;
  color?: string;
  style?: React.CSSProperties;
}> = ({ children, color = theme.accent, style }) => (
  <span
    style={{
      display: "inline-block",
      padding: "8px 16px",
      borderRadius: 999,
      border: `1px solid ${color}55`,
      background: `${color}14`,
      color,
      fontFamily: theme.fontMono,
      fontSize: 24,
      fontWeight: 600,
      ...style,
    }}
  >
    {children}
  </span>
);
