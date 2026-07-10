---
name: verify
description: Build/launch/drive recipe for verifying maglab changes end-to-end in a browser
---

# Verifying maglab

Static ES-module app — no build step. `npm test` runs the physics assertions,
but that is CI, not verification; verify by driving the app in a browser.

## Launch

```bash
python3 -m http.server 8123 &          # serve the repo root
curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8123/index.html   # expect 200
```

## Drive (headless Chromium via playwright-core)

Install `playwright-core` in a scratch dir (not the repo). Launch with the
pre-installed system browser — do NOT `playwright install`:

```js
const browser = await chromium.launch({
  executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',  // ls /opt/pw-browsers for the current version dir
  headless: true,
});
```

## Observing state

App state is not exposed on `window`; observe through the DOM and canvas:

- `#probeReadout`, `#forceReadout`, `#partReadout` — live data tiles (left panel)
- `#inspector` — shows "Select an object" when nothing is selected; param rows
  are `label.row` elements, find inputs via the row's `<span>` label text
- `#objlist .obj.sel` — selected object row
- Canvas interactions: `page.mouse.click/down/move/up` at coords from
  `#view`'s bounding box; world origin is the canvas centre at default view
- Collect `pageerror` + console `error` events across the whole drive —
  a single module-load error breaks everything silently otherwise

## Useful flows

- Add sources via `[data-add="wire"]` etc.; `#clearAll` resets the scene
- Force tile updates via rAF — wait ~200ms after a change before reading it
- Launch a particle: fill `#pSpeed`, click `#launch`; `#rateRow` becomes
  visible while the sim runs
