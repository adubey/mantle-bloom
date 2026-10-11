export async function asJson<T>(resp: Response): Promise<T> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.json() as Promise<T>;
}

export async function asBlob(resp: Response): Promise<Blob> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.blob();
}

export interface ColumnarStatsHistory {
  format: "columnar-v1";
  length: number;
  columns: Record<string, unknown[] | ColumnarMap>;
}

interface ColumnarMap {
  keys: string[];
  columns: Record<string, unknown[]>;
  states: Array<0 | 1 | 2>;
}

export function expandColumnarStatsHistory<T = Record<string, unknown>>(payload: ColumnarStatsHistory): T[] {
  if (payload.format !== "columnar-v1" || !Number.isInteger(payload.length) || payload.length < 0) {
    throw new Error("Invalid columnar stats history header");
  }

  const rows = Array.from({ length: payload.length }, () => ({} as Record<string, unknown>));
  for (const [field, column] of Object.entries(payload.columns)) {
    if (Array.isArray(column)) {
      if (column.length !== payload.length) throw new Error(`Invalid stats history column length: ${field}`);
      column.forEach((value, index) => {
        rows[index][field] = value;
      });
      continue;
    }

    const mapColumn = column as ColumnarMap;
    if (
      !Array.isArray(mapColumn.keys) || !Array.isArray(mapColumn.states) || mapColumn.states.length !== payload.length
    ) {
      throw new Error(`Invalid stats history map column: ${field}`);
    }
    for (const key of mapColumn.keys) {
      if (!Array.isArray(mapColumn.columns[key]) || mapColumn.columns[key].length !== payload.length) {
        throw new Error(`Invalid stats history map values: ${field}.${key}`);
      }
    }

    for (let index = 0; index < payload.length; index++) {
      const state = mapColumn.states[index];
      if (state === 0) {
        rows[index][field] = null;
      } else if (state === 1) {
        rows[index][field] = {};
      } else if (state === 2) {
        const value: Record<string, unknown> = {};
        for (const key of mapColumn.keys) {
          const item = mapColumn.columns[key][index];
          if (item !== null) value[key] = item;
        }
        rows[index][field] = value;
      } else {
        throw new Error(`Invalid stats history map state: ${field}`);
      }
    }
  }
  return rows as unknown as T[];
}
