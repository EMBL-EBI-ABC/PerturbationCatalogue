import { spawn } from "node:child_process";
import { writeFile } from "node:fs/promises";
import { pathToFileURL } from "node:url";

const [input, output] = process.argv.slice(2);
if (!input || !output) throw new Error("Usage: node render-poster.mjs INPUT.html OUTPUT.pdf");

const port = 9333;
const chrome = spawn("/snap/bin/chromium", [
  "--headless",
  "--disable-gpu",
  "--no-sandbox",
  "--allow-file-access-from-files",
  `--remote-debugging-port=${port}`,
  "--user-data-dir=/tmp/perturbation-catalogue-poster",
  "about:blank",
], { stdio: "ignore" });

const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
let target;
for (let attempt = 0; attempt < 50 && !target; attempt++) {
  try {
    target = await fetch(`http://127.0.0.1:${port}/json/new?about:blank`, { method: "PUT" }).then(response => response.json());
  } catch {
    await pause(100);
  }
}
if (!target) throw new Error("Chromium did not expose its print endpoint");

const socket = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((resolve, reject) => {
  socket.onopen = resolve;
  socket.onerror = reject;
});

let nextId = 0;
const pending = new Map();
const eventWaiters = new Map();
socket.onmessage = ({ data }) => {
  const message = JSON.parse(data);
  if (message.id) {
    const { resolve, reject } = pending.get(message.id);
    pending.delete(message.id);
    return message.error ? reject(new Error(message.error.message)) : resolve(message.result);
  }
  eventWaiters.get(message.method)?.splice(0).forEach(resolve => resolve(message.params));
};

const send = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++nextId;
  pending.set(id, { resolve, reject });
  socket.send(JSON.stringify({ id, method, params }));
});
const waitFor = method => new Promise(resolve => {
  const waiters = eventWaiters.get(method) ?? [];
  waiters.push(resolve);
  eventWaiters.set(method, waiters);
});

try {
  await send("Page.enable");
  const loaded = waitFor("Page.loadEventFired");
  await send("Page.navigate", { url: pathToFileURL(input).href });
  await loaded;
  await send("Runtime.evaluate", { expression: "document.fonts.ready", awaitPromise: true });
  const { data } = await send("Page.printToPDF", {
    printBackground: true,
    preferCSSPageSize: true,
    marginTop: 0,
    marginRight: 0,
    marginBottom: 0,
    marginLeft: 0,
  });
  await writeFile(output, Buffer.from(data, "base64"));
} finally {
  socket.close();
  chrome.kill();
}
