import type { ReactNode } from "react";

// One Generate World option: a checkbox with a title and a short grey explanation under it.
// Shared by the Random tab's terrain choices so each new choice isn't a copied block.
export default function GenerateOptionCheckbox({
  checked,
  onChange,
  disabled,
  title,
  children,
}: {
  checked: boolean;
  onChange: (checked: boolean) => void;
  disabled: boolean;
  title: string;
  children: ReactNode;
}) {
  return (
    <label style={{ display: "flex", alignItems: "flex-start", gap: 8, marginBottom: 16 }}>
      <input
        type="checkbox"
        checked={checked}
        onChange={(e) => onChange(e.target.checked)}
        disabled={disabled}
        style={{ marginTop: 2 }}
      />
      <span>
        {title}
        <span style={{ display: "block", fontSize: 11, color: "#999", marginTop: 3 }}>{children}</span>
      </span>
    </label>
  );
}
