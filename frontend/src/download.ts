// Browser file downloads, shared by the File dialog's saves and the recording's save dialog.

// How long an object URL outlives its download click. The click only starts the download, and
// the browser reads the blob after it returns, so revoking the URL straight away can cancel it.
const REVOKE_DELAY_MS = 60_000;

/** Start a browser download of `blob` as `filename`. */
export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.style.display = "none";
  document.body.appendChild(a);
  try {
    a.click();
  } finally {
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), REVOKE_DELAY_MS);
  }
}

/** Decode a base64 payload into a Blob without leaving the page. Fetching a `data:` URL would
 * be simpler, but it is a connection, and the app's `connect-src 'self'` policy refuses it, so
 * every recording's Save silently did nothing (#961). */
export function base64ToBlob(base64: string, mime: string): Blob {
  const binary = atob(base64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return new Blob([bytes], { type: mime });
}

export interface RecordedVideo {
  videoBase64: string;
  mime: string;
  seed: number;
  elapsed_years: number;
}

/** Download a finished recording's MP4. Throws if it can't, so the caller can keep the video. */
export function saveRecordedVideo(video: RecordedVideo): void {
  const blob = base64ToBlob(video.videoBase64, video.mime);
  downloadBlob(blob, `mantle-bloom-seed${video.seed}-animation-${Math.round(video.elapsed_years)}y.mp4`);
}
