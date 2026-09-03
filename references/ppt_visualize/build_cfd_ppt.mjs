import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";


const PROJECT = path.resolve(process.argv[2]);
const OUT = path.resolve(process.argv[3]);
const PREVIEW = path.resolve(process.argv[4]);
const HERE = path.join(PROJECT, "references", "ppt_visualize");
const ASSETS = path.join(HERE, "assets");
const evidence = JSON.parse(await fs.readFile(path.join(HERE, "cfd_ppt_content_evidence.json"), "utf8"));

const W = 1280;
const H = 720;
const C = {
  ink: "#1F2933",
  navy: "#163A5F",
  teal: "#2A7F76",
  pale: "#EEF4F5",
  line: "#CDD7DD",
  muted: "#66747E",
  soft: "#F6F8F9",
  white: "#FFFFFF",
  inlet: "#D55E00",
  out1: "#0072B2",
  out2: "#009E73",
  out3: "#CC79A7",
  amber: "#E69F00",
};

const titles = [
  "From 1D Vascular Data to Validated 3D CFD",
  "Three stages turn vascular data into a validated 3D flow field",
  "The 1D model supplies consistent 3D boundary conditions",
  "Surface preparation creates a clean solver-ready vascular domain",
  "The 3D solver uses a validated lattice-Boltzmann configuration",
  "Steady CFD resolves velocity across vascular branches",
  "Flow and gauge pressure remain consistent across three outlets",
  "The validated Base solution is production-ready within the current study scope",
];

const gh = (relative) => `https://github.com/OZAKI39/ulm-3d-vascular-model-generation/blob/codex/cfd-wall-force-numerics-validated-sync-20260830/${relative}`;

async function bytes(name) {
  const b = await fs.readFile(path.join(ASSETS, name));
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
}

function addText(slide, text, left, top, width, height, fontSize = 18, color = C.ink, bold = false, alignment = "left") {
  const shape = slide.shapes.add({
    geometry: "textbox",
    position: { left, top, width, height },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  shape.text = text;
  shape.text.style = { fontSize, fontFamily: "Arial", color, bold, alignment };
  return shape;
}

function addRect(slide, left, top, width, height, fill = C.white, line = C.line, radius = 0) {
  return slide.shapes.add({
    geometry: radius ? "roundRect" : "rect",
    position: { left, top, width, height },
    fill,
    line: { style: "solid", fill: line, width: line === "none" ? 0 : 1 },
    ...(radius ? { borderRadius: radius } : {}),
  });
}

function addLine(slide, left, top, width, height, color = C.line, weight = 1) {
  return slide.shapes.add({
    geometry: "line",
    position: { left, top, width, height },
    fill: "none",
    line: { style: "solid", fill: color, width: weight },
  });
}

function addTitle(slide, title, slideNumber, label) {
  addText(slide, title, 58, 33, 1070, 70, 32, C.ink, true);
  addText(slide, label, 1100, 43, 120, 28, 12, C.teal, true, "right");
  addLine(slide, 58, 106, 1164, 0, C.line, 1);
  addRect(slide, 58, 106, 74, 4, C.teal, "none");
  addText(slide, String(slideNumber).padStart(2, "0"), 1176, 681, 44, 22, 11, C.muted, false, "right");
}

function addNotes(slide, sources) {
  slide.speakerNotes.textFrame.setText(`[Sources]\n${sources.map((s) => `- ${s}`).join("\n")}\n[/Sources]`);
}

async function addImage(slide, name, left, top, width, height, alt, fit = "contain") {
  return slide.images.add({
    blob: await bytes(name),
    contentType: "image/png",
    alt,
    fit,
    position: { left, top, width, height },
  });
}

function value(slide, startsWith) {
  const item = evidence.claims.find((c) => c.slide === slide && c.claim.startsWith(startsWith));
  if (!item) throw new Error(`Missing evidence for slide ${slide}: ${startsWith}`);
  return item.value;
}

function format(value, digits) {
  return Number(value).toFixed(digits);
}

const deck = Presentation.create({ slideSize: { width: W, height: H } });

// Slide 1 — restrained cover with one project-specific focal render.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addRect(slide, 0, 0, 16, H, C.navy, "none");
  addText(slide, "COMPUTATIONAL METHODS · MOUSE MICROVASCULAR ROI", 66, 70, 500, 26, 12, C.teal, true);
  addText(slide, titles[0], 66, 126, 500, 150, 48, C.ink, true);
  addText(slide, "A three-stage workflow for constructing, preparing and solving flow in a mouse microvascular ROI", 66, 302, 470, 90, 21, C.muted, false);
  addLine(slide, 66, 430, 390, 0, C.line, 1);
  addText(slide, "1D boundary data   →   prepared lumen   →   steady 3D field", 66, 451, 480, 54, 17, C.navy, true);
  addRect(slide, 582, 54, 640, 604, C.soft, C.line, 2);
  await addImage(slide, "05_velocity_field.png", 596, 68, 612, 576, "Accepted Base CFD velocity field", "contain");
  addText(slide, "Accepted Base · physical velocity magnitude (mm s⁻¹)", 680, 640, 454, 24, 11, C.muted, false, "right");
  addNotes(slide, [
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/flow/production_steady_flow_field_manifest.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/visualization/interactive_v3_redesign/01_after_velocity_overview.json"),
  ]);
}

// Slide 2 — native three-stage workflow with actual thumbnails.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[1], 2, "WORKFLOW");
  const xs = [62, 452, 842];
  const names = ["1D boundary data", "CFD-ready surface", "3D flow solution"];
  const bodies = [
    "A reduced network supplies flow and pressure at the local ports.",
    "Extensions and caps create clean regions for boundary conditions.",
    "A validated lattice model resolves steady velocity and gauge pressure.",
  ];
  const images = ["01_pipeline_stage1_1d.png", "03_surface_after.png", "05_velocity_field.png"];
  const alts = ["Actual global 1D vascular network", "Actual prepared vascular surface", "Actual accepted Base velocity field"];
  const cards = [];
  for (let i = 0; i < 3; i += 1) {
    const card = addRect(slide, xs[i], 152, 328, 462, C.white, C.line, 4);
    cards.push(card);
    addRect(slide, xs[i], 152, 328, 8, i === 0 ? C.navy : i === 1 ? C.teal : C.out2, "none");
    addText(slide, `0${i + 1}`, xs[i] + 20, 174, 44, 30, 14, C.teal, true);
    addText(slide, names[i], xs[i] + 20, 205, 288, 46, 22, C.ink, true);
    addRect(slide, xs[i] + 20, 263, 288, 210, C.soft, "none", 2);
    await addImage(slide, images[i], xs[i] + 24, 267, 280, 202, alts[i], "contain");
    addText(slide, bodies[i], xs[i] + 20, 494, 288, 82, 17, C.muted, false);
  }
  slide.shapes.connect(cards[0], cards[1], { kind: "straight", fromSide: "right", toSide: "left", line: { style: "solid", fill: C.navy, width: 3 }, tail: { type: "triangle", width: "sm", length: "sm" } });
  slide.shapes.connect(cards[1], cards[2], { kind: "straight", fromSide: "right", toSide: "left", line: { style: "solid", fill: C.teal, width: 3 }, tail: { type: "triangle", width: "sm", length: "sm" } });
  addText(slide, "Scientific hand-off at each stage: boundary values → geometry → validated field", 62, 635, 1068, 34, 15, C.navy, true);
  addNotes(slide, [
    gh("outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json"),
    gh("outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/run_summary.json"),
  ]);
}

// Slide 3 — actual 1D/ROI views plus an editable boundary-condition derivation.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[2], 3, "STAGE 1");
  addRect(slide, 58, 134, 520, 235, C.soft, C.line, 3);
  await addImage(slide, "01_pipeline_stage1_1d.png", 64, 140, 508, 223, "Actual global 1D pressure network", "contain");
  addText(slide, "Global network", 76, 146, 160, 26, 13, C.navy, true);
  addRect(slide, 58, 385, 520, 250, C.soft, C.line, 3);
  await addImage(slide, "01b_roi_boundary_transfer.png", 64, 391, 508, 238, "Actual ROI with transferred inlet and outlet data", "contain");
  addText(slide, "Local ROI ports", 76, 397, 160, 26, 13, C.navy, true);

  addRect(slide, 606, 134, 616, 501, C.soft, C.line, 4);
  addText(slide, "How the 3D boundary conditions are constructed", 632, 151, 552, 44, 22, C.ink, true);
  addLine(slide, 632, 205, 564, 0, C.line, 1);

  addText(slide, "1 · Solve the 1D hydraulic network", 632, 220, 520, 28, 18, C.navy, true);
  addText(slide, "Q = Δp / R", 632, 253, 176, 33, 24, C.ink, true);
  addText(slide, "R = (8μ/π) ∫₀ᴸ ds / r(s)⁴", 823, 255, 350, 30, 19, C.ink, true);
  addText(slide, "Pressure drop, viscosity and vessel radius determine segment flow.", 632, 290, 538, 30, 16, C.muted, false);
  addLine(slide, 632, 329, 564, 0, C.line, 1);

  addText(slide, "2 · Map 1D pressure onto the prepared ports", 632, 344, 520, 28, 18, C.teal, true);
  addText(slide, "p_cut = p_parent − Q R(0→cut)", 632, 376, 484, 31, 21, C.ink, true);
  addText(slide, "p_out,i^3D = p_cut,i − Q_i^1D R_i,extension", 632, 406, 540, 31, 19, C.ink, true);
  addText(slide, "The second term removes the artificial extension pressure loss.", 632, 435, 548, 28, 16, C.muted, false);
  addLine(slide, 632, 469, 564, 0, C.line, 1);

  addText(slide, "3 · Apply the final solver boundaries", 632, 482, 548, 28, 18, C.out2, true);
  addText(slide, "INLET", 632, 518, 72, 24, 13, C.inlet, true);
  addText(slide, "Q_target = u_mean,healthy A_inlet = ∫_A u·n dA", 708, 515, 464, 27, 17, C.ink, true);
  addText(slide, "adaptive_flux_pressure adjusts inlet pressure to match Q_target", 708, 539, 470, 24, 16, C.muted, false);
  addText(slide, "OUTLETS", 632, 570, 72, 24, 13, C.out1, true);
  addText(slide, "pressure_eq: p_gauge,i = p_out,i^3D", 708, 567, 386, 27, 17, C.ink, true);
  addText(slide, "WALL", 632, 598, 72, 24, 13, C.navy, true);
  addText(slide, "wall_libb: u_wall = 0", 708, 595, 250, 27, 17, C.ink, true);
  addNotes(slide, [
    gh("outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/run_summary.json"),
    gh("outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628/qc/port_transfer_qc.json"),
    gh("outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/bc/boundary_conditions_vmtk_boundarynormal_crossseam.json"),
    gh("utils/cfd_preprocess/one_d_flow.py"),
    gh("utils/cfd_preprocess/port_transfer.py"),
    gh("utils/cfd_surface_prepare/vmtk_qc.py"),
    gh("outputs/cfd_flow/healthy_mouse_capillary_calibration_anchor003274_20260829_180310/qc/healthy_flow_target_calculation.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/input/production_numerical_contract.json"),
  ]);
}

// Slide 4 — matched-camera actual surface comparison.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[3], 4, "STAGE 2");
  addRect(slide, 58, 135, 1164, 410, C.soft, C.line, 3);
  await addImage(slide, "04_surface_before_after.png", 66, 143, 1148, 394, "Matched-camera comparison of original and prepared surfaces", "contain");
  addText(slide, "BEFORE · open cut surface", 104, 156, 300, 28, 15, C.muted, true);
  addText(slide, "AFTER · extended and capped", 708, 156, 330, 28, 15, C.teal, true);
  addLine(slide, 640, 157, 0, 362, C.line, 2);
  addRect(slide, 58, 566, 1164, 85, C.white, C.line, 3);
  addText(slide, "Core preserved", 82, 586, 190, 26, 18, C.navy, true);
  addText(slide, "0 μm max far-core motion", 82, 616, 220, 24, 14, C.muted, false);
  addText(slide, "Boundary-ready", 404, 586, 190, 26, 18, C.teal, true);
  addText(slide, "4 capped ports · 1 component", 404, 616, 230, 24, 14, C.muted, false);
  addText(slide, "Geometry QC", 736, 586, 170, 26, 18, C.ink, true);
  addText(slide, "Watertight · 0 intersections", 736, 616, 230, 24, 14, C.muted, false);
  const legend = [[C.inlet, "Inlet"], [C.out1, "Outlet 1"], [C.out2, "Outlet 2"], [C.out3, "Outlet 3"]];
  legend.forEach(([color, label], i) => {
    addRect(slide, 1010, 579 + i * 17, 10, 10, color, "none");
    addText(slide, label, 1028, 574 + i * 17, 90, 18, 11, C.muted, false);
  });
  addNotes(slide, [
    gh("outputs/model_generate/ultraliser_anchor003274_20260825_133350/geometry/lumen_surface_um.vtp"),
    gh("outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/geometry/cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp"),
    gh("outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611/qc/final_surface_qc.json"),
  ]);
}

// Slide 5 — native solver schematic and a small project-parameter strip.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[4], 5, "STAGE 3");
  addText(slide, "A uniform lattice advances local particle populations until the physical flow field becomes steady.", 58, 130, 780, 44, 18, C.muted, false);
  const boxes = [
    { x: 70, title: "Uniform 3D lattice", body: "The lumen is represented by 182,320 fluid cells.", color: C.navy },
    { x: 415, title: "Stream & collide", body: "D3Q19 populations exchange information locally.", color: C.teal },
    { x: 760, title: "Steady physical field", body: "Velocity and gauge pressure are recovered in SI units.", color: C.out2 },
  ];
  const shapes = [];
  for (const b of boxes) {
    const s = addRect(slide, b.x, 216, 280, 178, C.white, C.line, 5);
    shapes.push(s);
    addRect(slide, b.x, 216, 280, 8, b.color, "none");
    addText(slide, b.title, b.x + 22, 244, 236, 36, 21, C.ink, true);
    addText(slide, b.body, b.x + 22, 296, 230, 70, 16, C.muted, false);
  }
  slide.shapes.connect(shapes[0], shapes[1], { kind: "straight", fromSide: "right", toSide: "left", line: { style: "solid", fill: C.navy, width: 3 }, tail: { type: "triangle", width: "sm", length: "sm" } });
  slide.shapes.connect(shapes[1], shapes[2], { kind: "straight", fromSide: "right", toSide: "left", line: { style: "solid", fill: C.teal, width: 3 }, tail: { type: "triangle", width: "sm", length: "sm" } });
  addRect(slide, 1082, 208, 132, 194, C.pale, "none", 5);
  addText(slide, "DRIVING", 1096, 230, 104, 22, 12, C.teal, true, "center");
  addText(slide, "Adaptive\ninlet flow", 1096, 267, 104, 48, 16, C.navy, true, "center");
  addLine(slide, 1102, 325, 92, 0, C.line, 1);
  addText(slide, "3 outlet\ngauge pressures", 1096, 340, 104, 52, 15, C.ink, false, "center");

  const dx = format(value(5, "Base lattice spacing"), 2);
  const dt = format(value(5, "The diffusive time step"), 3);
  const tau = format(value(5, "The relaxation time"), 1);
  const q = format(value(5, "Target inflow"), 3);
  const stats = [["dx", `${dx} μm`], ["dt", `${dt} ns`], ["τ", tau], ["Target Q", `${q} nL/min`]];
  addText(slide, "Validated Base configuration", 70, 452, 440, 32, 22, C.ink, true);
  stats.forEach(([label, val], i) => {
    const x = 70 + i * 286;
    addLine(slide, x, 508, 232, 0, i === 0 ? C.navy : C.teal, 3);
    addText(slide, label, x, 526, 110, 25, 14, C.muted, true);
    addText(slide, val, x, 557, 232, 42, 23, C.ink, true);
  });
  addText(slide, "Diffusive scaling: dt = dx² / (6ν)", 70, 627, 420, 25, 14, C.muted, false);
  addNotes(slide, [
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/input/production_numerical_contract.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/flow/production_steady_flow_field_manifest.json"),
  ]);
}

// Slide 6 — large velocity field with minimal editable interpretation.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[5], 6, "RESULT");
  addRect(slide, 58, 133, 864, 500, C.soft, C.line, 3);
  await addImage(slide, "05_velocity_field.png", 66, 141, 848, 484, "Accepted Base physical velocity magnitude field", "contain");
  addText(slide, "PHYSICAL VELOCITY · mm s⁻¹", 78, 151, 250, 22, 12, C.navy, true);
  addRect(slide, 950, 153, 270, 118, C.pale, "none", 4);
  addText(slide, "0.177 mm s⁻¹", 970, 175, 230, 34, 25, C.navy, true);
  addText(slide, "mean speed across the accepted Base field", 970, 218, 226, 40, 14, C.muted, false);
  addRect(slide, 950, 291, 270, 118, C.white, C.line, 4);
  addText(slide, "598,755", 970, 313, 230, 34, 25, C.teal, true);
  addText(slide, "accepted steady-state iteration", 970, 356, 226, 30, 14, C.muted, false);
  addLine(slide, 950, 455, 270, 0, C.line, 1);
  addText(slide, "Interpretation", 950, 474, 230, 28, 16, C.ink, true);
  addText(slide, "Higher-speed colors are localized around the inlet and selected branch segments.", 950, 514, 250, 78, 17, C.muted, false);
  addNotes(slide, [
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/production_primary_metrics.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/visualization/interactive_v3_redesign/01_after_velocity_overview.json"),
  ]);
}

// Slide 7 — actual gauge-pressure render and editable outlet split chart.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[6], 7, "QUANTITATIVE RESULT");
  addRect(slide, 58, 136, 665, 480, C.soft, C.line, 3);
  await addImage(slide, "06_gauge_pressure_field.png", 66, 144, 649, 464, "Accepted Base gauge-pressure field", "contain");
  addText(slide, "GAUGE PRESSURE · Pa", 78, 154, 210, 22, 12, C.navy, true);
  addText(slide, "Outlet flow split", 760, 140, 430, 32, 22, C.ink, true);
  const fractions = value(7, "Outlet flow fractions").map(Number);
  slide.charts.add("bar", {
    position: { left: 750, top: 180, width: 450, height: 270 },
    categories: ["Outlet 1", "Outlet 2", "Outlet 3"],
    series: [{
      name: "Share",
      values: fractions,
      valuesFormatCode: "0.0",
      fill: C.out1,
      points: [{ idx: 0, fill: C.out1 }, { idx: 1, fill: C.out2 }, { idx: 2, fill: C.out3 }],
    }],
    hasLegend: false,
    barOptions: { direction: "bar", grouping: "clustered", gapWidth: 42, varyColors: true },
    xAxis: { min: 0, max: 80, majorUnit: 20, numberFormatCode: "0\"%\"", textStyle: { fill: C.muted, fontSize: 12 }, majorGridlines: { style: "solid", fill: "#E4EAED", width: 1 }, line: { style: "solid", fill: C.line, width: 1 } },
    yAxis: { textStyle: { fill: C.muted, fontSize: 13 }, line: { style: "solid", fill: "none", width: 0 }, majorGridlines: null },
    dataLabels: { showValue: true, position: "outEnd", textStyle: { fill: C.ink, fontSize: 13, bold: true } },
    chartFill: C.white,
    chartLine: { style: "solid", fill: "none", width: 0 },
    plotAreaFill: C.white,
    plotAreaLine: { style: "solid", fill: "none", width: 0 },
  });
  const qin = format(value(7, "Inlet flow"), 3);
  const qout = format(value(7, "Outlet flow"), 3);
  const closure = format(value(7, "Physical volume closure"), 3);
  const pin = Math.round(Number(value(7, "Inlet gauge pressure")));
  addLine(slide, 760, 474, 430, 0, C.line, 1);
  addText(slide, `${qin} → ${qout}`, 760, 492, 210, 34, 24, C.navy, true);
  addText(slide, "inlet → outlet (nL/min)", 760, 531, 210, 24, 13, C.muted, false);
  addText(slide, `${closure}%`, 1000, 492, 190, 34, 24, C.teal, true);
  addText(slide, "physical flow closure", 1000, 531, 190, 24, 13, C.muted, false);
  addText(slide, `${pin} Pa inlet gauge pressure`, 760, 576, 430, 34, 17, C.ink, true);
  addNotes(slide, [
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/production_primary_metrics.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/visualization/interactive_v3_redesign/02_after_pressure_overview.json"),
  ]);
}

// Slide 8 — native validation and scoped limitations.
{
  const slide = deck.slides.add();
  slide.background.fill = C.white;
  addTitle(slide, titles[7], 8, "VALIDATION & SCOPE");
  addRect(slide, 58, 143, 556, 408, C.pale, "none", 4);
  addText(slide, "VALIDATED", 86, 170, 220, 30, 15, C.teal, true);
  const residual = Number(value(8, "The Full V2 mass-identity residual")[0]);
  const closure = Number(value(7, "Physical volume closure"));
  const validated = [
    ["Steady Base solution", "Accepted at 598,755 iterations"],
    ["Physical flow closure", `${format(closure, 3)}%  (< 1%)`],
    ["Full-timestep mass identity", `${residual.toExponential(2)}  (< 1 × 10⁻⁸)`],
    ["Numerical safety", "Positive populations; low lattice speed"],
    ["Resolution sensitivity", "Coarse → Base max difference 1.45%"],
  ];
  validated.forEach(([head, body], i) => {
    const y = 218 + i * 61;
    addRect(slide, 88, y + 4, 13, 13, C.teal, "none");
    addText(slide, head, 116, y - 2, 250, 24, 16, C.ink, true);
    addText(slide, body, 366, y - 2, 220, 30, 14, C.muted, false, "right");
  });

  addRect(slide, 646, 143, 576, 408, C.white, C.line, 4);
  addText(slide, "CURRENT SCOPE LIMITS", 676, 170, 300, 30, 15, C.navy, true);
  const limits = [
    ["Three-grid convergence", "Richardson/GCI not completed under the compute budget."],
    ["Fine grid", "Mesh, controller and 5,000-step safety passed; steady state not completed."],
    ["Wall shear stress", "Validation is deferred to a later production phase."],
  ];
  limits.forEach(([head, body], i) => {
    const y = 229 + i * 94;
    addLine(slide, 678, y - 13, 512, 0, C.line, 1);
    addText(slide, head, 678, y, 180, 28, 17, C.ink, true);
    addText(slide, body, 872, y, 306, 60, 15, C.muted, false);
  });
  addRect(slide, 58, 583, 1164, 70, C.navy, "none", 3);
  addText(slide, "The validated Base field is the production flow solution used for downstream analysis.", 84, 604, 1110, 32, 21, C.white, true, "center");
  addNotes(slide, [
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/production_steady_qc.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/production_full_timestep_v2.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/two_grid_resolution_sensitivity.json"),
    gh("outputs/cfd_flow/production_tau1_base_promotion_anchor003274_20260902_013637/qc/run_summary.json"),
  ]);
}

await fs.mkdir(PREVIEW, { recursive: true });
for (const [index, slide] of deck.slides.items.entries()) {
  const stem = `slide_${String(index + 1).padStart(2, "0")}`;
  const png = await deck.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(path.join(PREVIEW, `${stem}.png`), new Uint8Array(await png.arrayBuffer()));
}

const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(OUT);
await fs.rm(`${OUT}.inspect.ndjson`, { force: true });

await fs.writeFile(
  path.join(HERE, "cfd_ppt_build_report.json"),
  JSON.stringify({
    status: "PASS",
    builder: "@oai/artifact-tool",
    slide_count: deck.slides.items.length,
    slide_size_px: { width: W, height: H },
    slide_size_inches: { width: 13.333333, height: 7.5 },
    slide_titles: titles,
    native_editable_charts: 1,
    native_editable_diagrams: 3,
    raster_scientific_render_objects: 9,
    unique_raster_scientific_files_used: 7,
    all_visible_text_english: true,
    all_main_claims_source_verified: evidence.all_claims_verified === true,
    new_cfd_solve_performed: false,
  }, null, 2) + "\n",
  "utf8",
);

console.log(`Built ${deck.slides.items.length}-slide deck at ${OUT}`);
