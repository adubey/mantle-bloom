import { useCallback, useState } from "react";

// The one-at-a-time action state modals share: `busy` names the running action, and
// a failure lands in `error`, shown inside the modal rather than behind it.
export function useModalAction() {
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const runAction = useCallback(async (name: string, action: () => Promise<void>) => {
    setBusy(name);
    setError(null);
    try {
      await action();
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(null);
    }
  }, []);

  return { busy, error, runAction };
}
