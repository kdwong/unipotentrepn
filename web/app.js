/*
 * Theta Atlas API contract
 * ------------------------
 * POST /api/calculate
 * request:  { "group": "so" | "mp", "orbit": [8, 6, 4, 4, 2, 2, 2, 2] }
 * response: {
 *   "orbit": [8, 6, 4, 4, 2, 2, 2, 2],
 *   "totals": [1, 2, 5, 6, 9, 12, 17, 22, 31],
 *   "tower": ["Mp(0)", "O(1)", "Mp(2)", ...],
 *   "final_form": {"p":15,"q":16},
 *   "allowed_bipartition_shapes": [{"primitive_pairs":[],"p":[3,2,1,1],"q":[4,2,1,1]}],
 *   "results": [{
 *     "k": 4,
 *     "label": "(1111000|00000000)",
 *     "concrete_path_count": 8,
 *     "groups": [{
 *       "pbp": {
 *         "p": ["***", "**", "c", "c"],
 *         "q": ["***d", "**", "s", "d"],
 *         "gamma": "B+",
 *         "primitive_pairs": [[2,3],[4,5],[8,9]],
 *         "extended_drc": [["***","**","c","c"], ["***da","**","s","d"]]
 *       },
 *       "outer_epsilon_labels": ["trivial"],
 *       "paths": [{
 *         "id": 1,
 *         "steps": ["O(0,1)[tt]", "Mp(2)(-1/2)", "..."],
 *         "twist_histories": [["tt", "tt", "dt", "dt", "tt"]],
 *         "outer_epsilon_label": "trivial"
 *       }]
 *     }]
 *   }]
 * }
 *
 * Each string in p and q is one COLUMN of the painted diagram. Characters run
 * top-to-bottom inside a column; columns are rendered left-to-right, top-aligned.
 * The normalizer below also accepts common snake/camel-case aliases.
 */

(function () {
  "use strict";

  const configuredApiUrl =
    document.querySelector('meta[name="theta-api-url"]')?.content.trim() || "";
  const API_URL = configuredApiUrl || "/api/calculate";
  const form = document.querySelector("#orbit-form");
  const groupInputs = [...document.querySelectorAll('input[name="group"]')];
  const programGroup = document.querySelector("#program-group");
  const orbitInput = document.querySelector("#orbit-input");
  const orbitInputWrap = document.querySelector(".orbit-input-wrap");
  const orbitError = document.querySelector("#orbit-error");
  const calculateButton = document.querySelector("#calculate-button");
  const buttonLabel = calculateButton.querySelector(".button-label");
  const statusRegion = document.querySelector("#status-region");
  const loadingView = document.querySelector("#loading-view");
  const resultsView = document.querySelector("#results-view");
  const resultsTitle = document.querySelector("#results-title");
  const resultsSubtitle = document.querySelector("#results-subtitle");
  const resultStats = document.querySelector("#result-stats");
  const towerElement = document.querySelector("#tower");
  const finalFormElement = document.querySelector("#final-form");
  const kNavigationLabel = document.querySelector("#k-navigation-label");
  const kLinks = document.querySelector("#k-links");
  const kResults = document.querySelector("#k-results");
  const pbpTemplate = document.querySelector("#pbp-group-template");
  const runtimeLabel = document.querySelector(".runtime-label");
  const DIRECT_FILE_MESSAGE =
    "This HTML file cannot run the calculation by itself. Close this tab, then double-click start_theta2_pbp_web.bat in the BMSZ folder.";

  let activeRequest = null;

  if (configuredApiUrl) {
    runtimeLabel.textContent = "Online computation";
  }

  form.addEventListener("submit", handleSubmit);
  groupInputs.forEach((input) => input.addEventListener("change", handleGroupChange));

  document.querySelectorAll(".example-chip").forEach((button) => {
    button.addEventListener("click", () => {
      orbitInput.value = button.dataset.orbit;
      clearError();
      form.requestSubmit();
    });
  });

  orbitInput.addEventListener("input", clearError);
  updateGroupCopy();

  if (window.location.protocol === "file:") {
    calculateButton.disabled = true;
    buttonLabel.textContent = "Start server first";
    showError(DIRECT_FILE_MESSAGE);
  }

  async function handleSubmit(event) {
    event.preventDefault();

    if (window.location.protocol === "file:") {
      showError(DIRECT_FILE_MESSAGE);
      return;
    }

    let orbit;
    try {
      orbit = parseOrbit(orbitInput.value);
    } catch (error) {
      showError(error.message);
      return;
    }
    const group = selectedGroup();
    if (activeRequest) activeRequest.abort();
    const requestController = new AbortController();
    activeRequest = requestController;
    setLoading(true, orbit);

    try {
      const response = await fetch(API_URL, {
        method: "POST",
        headers: { "Content-Type": "application/json", Accept: "application/json" },
        body: JSON.stringify({ group, orbit }),
        signal: requestController.signal,
      });

      const payload = await readPayload(response);
      if (!response.ok) {
        const reason =
          payload?.error?.message ||
          (typeof payload?.error === "string" ? payload.error : "") ||
          payload?.message ||
          `Calculation failed (${response.status}).`;
        throw new Error(reason);
      }

      if (group === "mp" && payload?.group !== "mp") {
        throw new Error(
          "The online Mp calculator service is still updating. Please try again shortly.",
        );
      }

      const result = normalizeResponse(payload, orbit, group);
      renderResults(result);
      setLoading(false);
      resultsView.hidden = false;
      announce(`Calculation complete for ${result.groupLabel}. Found ${result.results.length} fine K-types.`);
      resultsView.scrollIntoView({ behavior: reducedMotion() ? "auto" : "smooth", block: "start" });
    } catch (error) {
      if (error.name === "AbortError") return;
      setLoading(false);
      showError(humanizeFetchError(error));
    } finally {
      if (activeRequest === requestController) activeRequest = null;
    }
  }

  function selectedGroup() {
    return groupInputs.find((input) => input.checked)?.value === "mp" ? "mp" : "so";
  }

  function handleGroupChange() {
    if (activeRequest) activeRequest.abort();
    activeRequest = null;
    setLoading(false);
    resultsView.hidden = true;
    clearError();
    updateGroupCopy();
    orbitInput.focus();
  }

  function updateGroupCopy() {
    const isMp = selectedGroup() === "mp";
    programGroup.textContent = isMp ? "G = Mp(2n, ℝ)" : "G = SO(n,n+1)";
  }

  function parseOrbit(value) {
    const cleaned = value
      .trim()
      .replace(/[()[\]{}]/g, " ")
      .replace(/O\s*(?:\^\s*)?(?:\{\s*)?(?:\\?vee|∨)(?:\s*\})?\s*=?/gi, " ")
      .trim();

    if (!cleaned) throw new Error("Enter at least one row size.");
    if (!/^[\d,;\s]+$/.test(cleaned)) {
      throw new Error("Use only even integers, commas, or spaces.");
    }

    const hasSeparator = /[\s,;]/.test(cleaned);
    const compactDigits = [...cleaned].map(Number);
    const looksLikeCompactOrbit =
      !hasSeparator &&
      compactDigits.length > 1 &&
      compactDigits.every((row) => row > 0 && row % 2 === 0) &&
      compactDigits.every((row, index) => index === 0 || compactDigits[index - 1] >= row);
    const orbit = looksLikeCompactOrbit
      ? compactDigits
      : cleaned.split(/[\s,;]+/).filter(Boolean).map(Number);

    if (!orbit.length || orbit.some((row) => !Number.isSafeInteger(row) || row <= 0)) {
      throw new Error("Every row size must be a positive integer.");
    }
    if (orbit.some((row) => row % 2 !== 0)) {
      throw new Error("Good parity here requires every row size to be even.");
    }
    if (orbit.some((row, index) => index > 0 && orbit[index - 1] < row)) {
      throw new Error("List the row sizes in decreasing order.");
    }
    return orbit;
  }

  async function readPayload(response) {
    const text = await response.text();
    if (!text) return {};
    try {
      return JSON.parse(text);
    } catch (_error) {
      if (!response.ok) return { error: text };
      throw new Error("The calculator returned an unreadable response.");
    }
  }

  function humanizeFetchError(error) {
    if (window.location.protocol === "file:") return DIRECT_FILE_MESSAGE;
    if (error instanceof TypeError) {
      if (configuredApiUrl) {
        return "The online calculation service could not be reached. Please try again in a moment.";
      }
      return "The local Python calculator is no longer running. Close this tab, then double-click start_theta2_pbp_web.bat in the BMSZ folder.";
    }
    return error.message || "The calculation could not be completed.";
  }

  function setLoading(isLoading, orbit = []) {
    calculateButton.disabled = isLoading;
    calculateButton.classList.toggle("is-loading", isLoading);
    buttonLabel.textContent = isLoading ? "Working" : "Calculate";
    loadingView.hidden = !isLoading;

    if (isLoading) {
      resultsView.hidden = true;
      clearError();
      announce(`Calculating the orbit (${orbit.join(", ")}). Larger orbits may take a moment.`);
    }
  }

  function showError(message) {
    orbitError.textContent = message;
    orbitInputWrap.classList.add("invalid");
    orbitInput.setAttribute("aria-invalid", "true");
    announce(message);
  }

  function clearError() {
    orbitError.textContent = "";
    orbitInputWrap.classList.remove("invalid");
    orbitInput.removeAttribute("aria-invalid");
  }

  function announce(message) {
    statusRegion.textContent = "";
    window.setTimeout(() => {
      statusRegion.textContent = message;
    }, 20);
  }

  function normalizeResponse(payload, requestedOrbit, requestedGroup) {
    if (!payload || typeof payload !== "object") {
      throw new Error("The calculator returned no result.");
    }

    const groupKind = normalizeGroupKind(payload.group ?? requestedGroup);
    const groupLabel = String(
      payload.group_label
        ?? payload.groupLabel
        ?? (groupKind === "mp" ? "Mp(2n, ℝ)" : "SO(n,n+1)"),
    );
    const orbit = toNumberArray(payload.orbit || payload.partition || requestedOrbit);
    const totals = toNumberArray(payload.totals || payload.chain_totals || payload.chainTotals || []);
    const finalForm = groupKind === "so"
      ? normalizeFinalForm(payload.final_form || payload.finalForm, orbit)
      : null;
    const ambientGroup = String(
      payload.ambient_group
        ?? payload.ambientGroup
        ?? `Sp(${orbit.reduce((sum, row) => sum + row, 0)}, ℂ)`,
    );
    const finalGroupLabel = String(
      payload.final_group_label
        ?? payload.finalGroupLabel
        ?? (groupKind === "mp"
          ? `Mp(${orbit.reduce((sum, row) => sum + row, 0)}, ℝ)`
          : `O(${finalForm.p}, ${finalForm.q})`),
    );
    const expectedBipartition = normalizeExpectedBipartition(
      payload.expected_bipartition || payload.expectedBipartition,
    );
    const allowedBipartitionShapes = normalizeAllowedShapes(
      payload.allowed_bipartition_shapes || payload.allowedBipartitionShapes,
      expectedBipartition,
    );
    const tower = normalizeTower(payload.tower, totals);
    const sourceResults = payload.results || payload.k_results || payload.kResults || payload.results_by_k;
    const results = normalizeKResults(sourceResults, finalForm, groupKind);

    if (!results.length) {
      throw new Error("The calculation returned no fine K-type results for this orbit.");
    }

    return {
      groupKind,
      groupLabel,
      orbit,
      totals,
      finalForm,
      finalGroupLabel,
      ambientGroup,
      expectedBipartition,
      allowedBipartitionShapes,
      tower,
      results,
    };
  }

  function normalizeGroupKind(value) {
    return String(value || "so").toLowerCase() === "mp" ? "mp" : "so";
  }

  function normalizeExpectedBipartition(value) {
    if (!value || typeof value !== "object") return { p: [], q: [] };
    return {
      p: toNumberArray(value.p || value.P || []),
      q: toNumberArray(value.q || value.Q || []),
    };
  }

  function normalizeAllowedShapes(value, fallback) {
    if (!Array.isArray(value) || !value.length) {
      return [{ primitivePairsLabel: "∅", p: fallback.p, q: fallback.q }];
    }
    return value.map((shape) => ({
      primitivePairsLabel: String(shape.primitive_pairs_label ?? shape.primitivePairsLabel ?? "∅"),
      p: toNumberArray(shape.p || shape.P || []),
      q: toNumberArray(shape.q || shape.Q || []),
    }));
  }

  function normalizeFinalForm(value, orbit) {
    const n = orbit.reduce((sum, row) => sum + row, 0) / 2;
    if (Array.isArray(value)) return { p: Number(value[0]), q: Number(value[1]) };
    if (value && typeof value === "object") {
      return { p: Number(value.p), q: Number(value.q) };
    }
    return { p: n, q: n + 1 };
  }

  function normalizeTower(value, totals) {
    if (Array.isArray(value) && value.length) {
      return value.map((stage, index) => {
        if (typeof stage === "string") return { label: stage };
        if (Array.isArray(stage)) return { label: `${stage[0]}(${stage[1]})` };
        const kind = stage.kind || stage.type || stage.group || (index % 2 ? "O" : "Mp");
        const size = stage.size ?? stage.total ?? stage.dimension ?? "";
        return { label: stage.label || `${kind}(${size})` };
      });
    }

    return [
      { label: "Mp(0)" },
      ...totals.map((total, index) => ({ label: `${index % 2 === 0 ? "O" : "Mp"}(${total})` })),
    ];
  }

  function normalizeKResults(value, finalForm, groupKind) {
    let entries;
    if (Array.isArray(value)) {
      entries = value.map((section) => [section?.label, section]);
    } else if (value && typeof value === "object") {
      entries = Object.entries(value);
    } else {
      entries = [];
    }

    return entries
      .map(([key, section], index) =>
        normalizeKSection(section || {}, key, index, finalForm, groupKind),
      )
      .sort((a, b) => a.k - b.k);
  }

  function normalizeKSection(section, key, index, finalForm, groupKind) {
    const label = String(section.label || section.k_type || section.kType || key || "");
    const explicitK = section.k ?? section.degree;
    const k = Number.isFinite(Number(explicitK)) ? Number(explicitK) : inferK(label, index);
    const groupsSource = section.groups || section.pbp_groups || section.pbpGroups || section.bipartitions || [];
    const groups = Array.isArray(groupsSource)
      ? groupsSource.map((group, groupIndex) => normalizeGroup(group, groupIndex))
      : Object.entries(groupsSource).map(([groupKey, group], groupIndex) =>
          normalizeGroup({ key: groupKey, ...(group || {}) }, groupIndex),
        );

    const concretePathCount = numberOr(
      section.concrete_path_count
        ?? section.concretePathCount
        ?? section.certified_count
        ?? section.certifiedCount
        ?? section.path_count
        ?? section.pathCount,
      uniquePathCount(groups),
    );
    const distinctPbpCount = numberOr(
      section.distinct_pbp_count
        ?? section.distinctPbpCount
        ?? section.distinct_so_pbp_count
        ?? section.distinctSoPbpCount,
      groups.length,
    );
    const oWedgeDegrees = toNumberArray(section.o_wedge_degrees || section.oWedgeDegrees || []);

    return {
      k,
      label: label || (groupKind === "so" ? constructLabel(k, finalForm) : `k = ${k}`),
      title: String(section.title ?? section.heading ?? ""),
      degreeDescription: String(
        section.degree_description ?? section.degreeDescription ?? "",
      ),
      concretePathCount,
      distinctPbpCount,
      oWedgeDegrees,
      middleDegree: Boolean(section.middle_degree ?? section.middleDegree),
      groups,
    };
  }

  function normalizeGroup(group, index) {
    const pbp = normalizePbp(group.pbp || group.painted_bipartition || group.paintedBipartition || group);
    const associatedCycle = normalizeAssociatedCycle(
      group.associated_cycle || group.associatedCycle || {},
    );
    const pathsSource = group.paths || group.chains || [];
    const paths = Array.isArray(pathsSource)
      ? pathsSource.map((path, pathIndex) => normalizePath(path, pathIndex))
      : Object.values(pathsSource).map((path, pathIndex) => normalizePath(path, pathIndex));
    const epsilonSource = group.outer_epsilon_labels ?? group.outerEpsilonLabels ?? [];
    const outerEpsilons = Array.isArray(epsilonSource)
      ? epsilonSource.map(String)
      : [String(epsilonSource)];
    return {
      id: group.id || group.key || `pbp-${index + 1}`,
      pbp,
      associatedCycle,
      paths,
      outerEpsilons,
    };
  }

  function normalizePbp(value) {
    const plain = value.plain_drc || value.plainDrc || value.drc;
    const p = stringArray(value.p || value.P || plain?.[0] || []);
    const q = stringArray(value.q || value.Q || plain?.[1] || []);
    const primitivePairsSource = value.primitive_pairs || value.primitivePairs || [];
    const primitivePairs = Array.isArray(primitivePairsSource)
      ? primitivePairsSource
          .filter((pair) => Array.isArray(pair) && pair.length >= 2)
          .map((pair) => [Number(pair[0]), Number(pair[1])])
      : [];

    const gammaValue = value.gamma ?? value.sign;
    return {
      p,
      q,
      gamma: gammaValue == null ? "" : String(gammaValue),
      parameterType: String(
        value.parameter_type ?? value.parameterType ?? value.type ?? "",
      ),
      primitivePairs,
      primitivePairsLabel: String(
        value.primitive_pairs_label
          ?? value.primitivePairsLabel
          ?? (primitivePairs.length
            ? `{${primitivePairs.map((pair) => `(${pair[0]},${pair[1]})`).join(", ")}}`
            : "∅"),
      ),
      raw: value.extended_drc || value.extendedDrc || value.raw_extended_drc || value.rawExtendedDrc || null,
      specialReference: value.special_reference_extended_drc || value.specialReferenceExtendedDrc || null,
    };
  }

  function normalizeAssociatedCycle(value) {
    const termsSource = Array.isArray(value.terms) ? value.terms : [];
    const terms = termsSource.map((term) => ({
      coefficient: numberOr(term.coefficient ?? term.multiplicity, 1),
      markedRows: stringArray(term.marked_rows || term.markedRows || []),
      underlyingPartition: toNumberArray(
        term.underlying_partition || term.underlyingPartition || [],
      ),
    }));
    return {
      termCount: numberOr(value.term_count ?? value.termCount, terms.length),
      totalMultiplicity: numberOr(
        value.total_multiplicity ?? value.totalMultiplicity,
        terms.reduce((total, term) => total + term.coefficient, 0),
      ),
      terms,
    };
  }

  function subsetLabel(indices) {
    return `{${indices.join(", ")}}`;
  }

  function normalizePath(value, index) {
    if (typeof value === "string") {
      const steps = splitChain(value);
      return {
        id: index + 1,
        number: index + 1,
        name: `Path ${index + 1}`,
        steps,
        histories: [],
        realizations: [{ history: [], historyText: "", steps }],
        outerEpsilon: "—",
      };
    }

    const rawSteps = value.steps || value.chain_steps || value.chainSteps || value.chain || [];
    const steps = Array.isArray(rawSteps) ? rawSteps.map(String) : splitChain(String(rawSteps));
    const historiesSource =
      value.twist_histories || value.twistHistories || value.histories || value.concrete_twists || [];
    const histories = normalizeHistories(historiesSource);
    const realizationSource =
      value.concrete_realizations || value.concreteRealizations || [];
    const realizations = normalizeRealizations(realizationSource, histories, steps);
    const number = numberOr(value.number ?? value.path_number ?? value.pathNumber, index + 1);
    const id = value.id ?? number;
    const subsetIndices = toNumberArray(
      value.subset_indices || value.subsetIndices || [],
    );
    const pathSubsetLabel = String(
      value.subset_label
        ?? value.subsetLabel
        ?? (subsetIndices.length ? subsetLabel(subsetIndices) : ""),
    );

    return {
      id,
      number,
      name: value.name || value.title || `Path ${pathSubsetLabel || number}`,
      subsetIndices,
      subsetLabel: pathSubsetLabel,
      targetLabel: String(
        value.target_label
          ?? value.targetLabel
          ?? value.final?.exact_label
          ?? value.final?.exactLabel
          ?? ""
      ),
      leftDegree: numberOr(
        value.left_degree
          ?? value.leftDegree
          ?? value.final?.left_degree
          ?? value.final?.leftDegree,
        Number.POSITIVE_INFINITY,
      ),
      steps,
      histories,
      realizations,
      outerEpsilon: String(
        value.outer_epsilon_label ?? value.outerEpsilonLabel ?? ""
      ),
    };
  }

  function normalizeHistories(value) {
    if (value == null || value === "") return [];
    if (!Array.isArray(value)) return [[String(value)]];
    if (!value.length) return [];
    if (value.every((item) => typeof item === "string")) {
      const strings = value.map(String);
      const lookLikeSeparateHistories = strings.some((item) => item.includes("->") || item.includes("→"));
      return lookLikeSeparateHistories ? strings.map((item) => [item]) : [strings];
    }
    return value.map((history) => (Array.isArray(history) ? history.map(String) : [String(history)]));
  }

  function normalizeRealizations(value, fallbackHistories, fallbackSteps) {
    if (Array.isArray(value) && value.length) {
      return value.map((realization, index) => {
        const historySource =
          realization?.twist_history
          ?? realization?.twistHistory
          ?? realization?.history
          ?? fallbackHistories[index]
          ?? [];
        const history = Array.isArray(historySource)
          ? historySource.map(String)
          : [String(historySource)];
        const historyText = String(
          realization?.twist_history_text
          ?? realization?.twistHistoryText
          ?? realization?.historyText
          ?? history.join(" → "),
        );
        const rawSteps = realization?.steps || realization?.chain_steps || fallbackSteps;
        const steps = Array.isArray(rawSteps)
          ? rawSteps.map(String)
          : splitChain(String(rawSteps || ""));
        const exactCycleSource =
          realization?.exact_associated_cycle
          ?? realization?.exactAssociatedCycle
          ?? null;
        return {
          history,
          historyText,
          steps,
          exactAssociatedCycle: exactCycleSource
            ? normalizeAssociatedCycle(exactCycleSource)
            : null,
          pbpCycleRelation: String(
            realization?.pbp_cycle_relation
            ?? realization?.pbpCycleRelation
            ?? "",
          ),
          pbpCycleTwist: toNumberArray(
            realization?.pbp_cycle_twist
            ?? realization?.pbpCycleTwist
            ?? [],
          ),
          pbpCycleRelationLabel: String(
            realization?.pbp_cycle_relation_label
            ?? realization?.pbpCycleRelationLabel
            ?? "",
          ),
        };
      });
    }
    if (fallbackHistories.length) {
      return fallbackHistories.map((history) => ({
        history,
        historyText: history.join(" → "),
        steps: fallbackSteps,
        exactAssociatedCycle: null,
        pbpCycleRelation: "",
        pbpCycleTwist: [],
        pbpCycleRelationLabel: "",
      }));
    }
    return [{
      history: [],
      historyText: "",
      steps: fallbackSteps,
      exactAssociatedCycle: null,
      pbpCycleRelation: "",
      pbpCycleTwist: [],
      pbpCycleRelationLabel: "",
    }];
  }

  function splitChain(value) {
    if (!value) return [];
    return value.split(/\s+→\s+|\s+->\s+/).map((step) => step.trim()).filter(Boolean);
  }

  function inferK(label, fallback) {
    const left = label.replace(/[()\s]/g, "").split("|")[0];
    if (/^1*0*$/.test(left)) return (left.match(/^1*/) || [""])[0].length;
    const exponent = label.match(/1\s*\^\s*\{?(\d+)/);
    return exponent ? Number(exponent[1]) : fallback;
  }

  function constructLabel(k, finalForm) {
    const leftSize = Math.max(0, Number(finalForm.p) || 0);
    const rightSize = Math.max(0, Number(finalForm.q) || 0);
    return `(${"1".repeat(k)}${"0".repeat(Math.max(0, leftSize - k))}|${"0".repeat(rightSize)})`;
  }

  function renderResults(data) {
    const orbitText = `(${data.orbit.join(", ")})`;
    const dualMark = document.createElement("sup");
    dualMark.textContent = "∨";
    resultsTitle.replaceChildren(
      document.createTextNode("O"),
      dualMark,
      document.createTextNode(` = ${orbitText}`),
    );
    const shapeNote = data.allowedBipartitionShapes.length === 1
      ? ` · ℘ = ∅ shape P = (${data.allowedBipartitionShapes[0].p.join(", ")}), Q = (${data.allowedBipartitionShapes[0].q.join(", ")})`
      : ` · ${data.allowedBipartitionShapes.length} legitimate ℘-dependent PBP shapes`;
    resultsSubtitle.textContent = `A good-parity orbit in ${data.ambientGroup}${shapeNote}`;

    const totalPaths = data.results.reduce((sum, section) => sum + section.concretePathCount, 0);
    const totalPbps = data.results.reduce((sum, section) => sum + section.distinctPbpCount, 0);
    resultStats.replaceChildren(
      makeStat(data.results.length, "fine K-types"),
      makeStat(totalPaths, "concrete theta paths"),
      makeStat(totalPbps, "painted bipartitions"),
    );

    renderTower(data.tower);
    finalFormElement.textContent = `final real group  ${data.finalGroupLabel}`;
    kNavigationLabel.textContent = data.groupKind === "mp"
      ? "Fine Mp K-type"
      : "Fine K-type";

    kLinks.replaceChildren(...data.results.map((section) => makeKLink(section, data.groupKind)));
    kResults.replaceChildren(
      ...data.results.map((section) => renderKSection(section, data.groupKind)),
    );
  }

  function makeStat(value, label) {
    const stat = element("div", "stat");
    stat.append(element("strong", "", String(value)), element("span", "", label));
    return stat;
  }

  function renderTower(tower) {
    const stages = tower.length ? tower : [{ label: "Mp(0)" }];
    towerElement.replaceChildren(
      ...stages.map((stage, index) => {
        const wrapper = element("div", "tower-stage");
        wrapper.setAttribute("role", "listitem");
        wrapper.append(element("div", "tower-node", stage.label));
        if (index < stages.length - 1) wrapper.append(element("span", "tower-arrow", "→"));
        return wrapper;
      }),
    );
  }

  function makeKLink(section, groupKind) {
    const link = element("a", "k-link");
    link.href = `#k-${section.k}`;
    link.setAttribute("aria-label", `Jump to fine K-type k equals ${section.k}`);
    if (groupKind === "mp") {
      link.textContent = String(section.k);
    } else {
      const one = document.createTextNode("1");
      const sup = element("sup", "", String(section.k));
      link.append(one, sup);
    }
    return link;
  }

  function renderKSection(section, groupKind) {
    const article = element("article", "k-result");
    article.id = `k-${section.k}`;

    const header = element("header", "k-result-header");
    const type = element("div", "k-type-label");
    type.append(element("span", "k-ordinal", `k = ${section.k}`));
    type.append(element("code", "exact-k-label", section.label));

    const counts = element("p", "k-counts");
    counts.append(
      strongText(section.distinctPbpCount),
      document.createTextNode(` painted bipartition${plural(section.distinctPbpCount)} · `),
      strongText(section.concretePathCount),
      document.createTextNode(` concrete theta path${plural(section.concretePathCount)}`),
    );
    counts.title = `Exact label ${section.label}`;
    header.append(type, counts);
    article.append(header);

    if (groupKind === "mp" && section.degreeDescription) {
      article.append(element("p", "multi-pbp-note", section.degreeDescription));
    } else if (section.oWedgeDegrees.length) {
      const degreeText = section.oWedgeDegrees.join(" and ");
      article.append(
        element(
          "p",
          "multi-pbp-note",
          section.middleDegree
            ? `The exact right-trivial left degree is ${degreeText}. Here p = 2k, so there is only one complementary degree; every twist history is nevertheless computed.`
            : `The exact right-trivial targets have left O(p)-degrees ${degreeText}. Paths to both target degrees are computed before their painted bipartitions are grouped; no count is multiplied.`,
        ),
      );
    }

    article.append(
      element(
        "p",
        "multi-pbp-note",
        groupKind === "mp"
          ? "The 2^r path subsets index the r rows of the dual orbit from bottom to top; the selected half-row lengths count the -1/2 entries and sum to k. These are the painted bipartitions reached by the fine-K-type paths, not an enumeration of every type-M extended PBP. The primitive-pair set ℘ uses the separate BMSZ top-to-bottom indexing, and repeated equal rows can make different paths share one painted bipartition."
          : "Path subsets index the rows of the dual orbit from bottom to top. The cycle shown with each painted bipartition is the canonical PBP reference; every expanded path shows the exact O(p,q)-extension cycle and whether it is unchanged or twisted by ⊗ (1,1).",
      ),
    );

    const groupedPathOccurrences = section.groups.reduce(
      (total, group) => total + group.paths.length,
      0,
    );
    if (groupedPathOccurrences > section.concretePathCount) {
      article.append(
        element(
          "p",
          "multi-pbp-note",
          "A concrete theta path can reach more than one painted bipartition. Its row-subset label is repeated in every group it reaches.",
        ),
      );
    }

    if (!section.groups.length) {
      article.append(element("div", "no-results", "No painted bipartitions occur for this K-type."));
      return article;
    }

    const groups = element("div", "pbp-groups");
    section.groups.forEach((group, index) =>
      groups.append(renderPbpGroup(group, index, section.k, groupKind)),
    );
    article.append(groups);
    return article;
  }

  function renderPbpGroup(group, index, k, groupKind) {
    const fragment = pbpTemplate.content.cloneNode(true);
    const card = fragment.querySelector(".pbp-group");
    const number = String(index + 1).padStart(2, "0");
    fragment.querySelector(".pbp-number").textContent = number;
    fragment.querySelector(".pbp-overline").textContent = `${group.paths.length} concrete theta path${plural(group.paths.length)} for this ${groupKind === "mp" ? "Mp" : "SO"} parameter`;
    fragment.querySelector(".pbp-title").textContent = `Painted bipartition ${number}`;

    const tags = fragment.querySelector(".pbp-tags");
    if (group.pbp.gamma) {
      tags.append(element("span", "tag gamma", `γ ${group.pbp.gamma}`));
    } else if (group.pbp.parameterType) {
      tags.append(element("span", "tag gamma", `type ${group.pbp.parameterType}`));
    }
    tags.append(element("span", "tag wp", `℘ ${group.pbp.primitivePairsLabel}`));
    if (groupKind === "so") {
      tags.append(
        element("span", "tag outer", `BMSZ outer ε · ${group.outerEpsilons.join(", ") || "none"}`),
      );
    }

    renderTableau(fragment.querySelector(".tableau-p"), group.pbp.p, "P");
    renderTableau(fragment.querySelector(".tableau-q"), group.pbp.q, "Q");
    const cyclePanel = fragment.querySelector(".associated-cycle-panel");
    if (groupKind === "so") {
      cyclePanel.setAttribute(
        "aria-label",
        "Painted-bipartition reference associated cycle",
      );
      cyclePanel.querySelector(".cycle-kicker").textContent =
        "Painted-bipartition cycle";
      cyclePanel.querySelector(".associated-cycle-heading h5").textContent =
        "Canonical marked diagrams";
      cyclePanel.querySelector(".associated-cycle-heading").after(
        element(
          "p",
          "cycle-reference-note",
          "Each exact O(p,q)-path cycle below is checked against this reference and is labelled either as displayed or ⊗ (1,1).",
        ),
      );
    }
    renderAssociatedCycle(
      cyclePanel,
      group.associatedCycle,
    );

    fragment.querySelector(".paths-heading h5").textContent =
      `${group.paths.length} concrete theta path${plural(group.paths.length)}`;
    const list = fragment.querySelector(".paths-list");
    [...group.paths]
      .sort((left, right) => left.leftDegree - right.leftDegree || left.number - right.number)
      .forEach((path, pathIndex) =>
        list.append(renderPath(path, pathIndex, k, index, groupKind)),
      );

    const toggle = fragment.querySelector(".toggle-paths");
    toggle.addEventListener("click", () => {
      const details = [...list.querySelectorAll("details")];
      const shouldOpen = details.some((detail) => !detail.open);
      details.forEach((detail) => {
        detail.open = shouldOpen;
      });
      updateToggleLabel(toggle, details);
    });
    list.addEventListener("toggle", () => updateToggleLabel(toggle, [...list.querySelectorAll("details")]), true);

    card.dataset.pbp = group.id;
    return fragment;
  }

  function renderAssociatedCycle(container, cycle) {
    const count = cycle.termCount;
    container.querySelector(".cycle-count").textContent =
      `${count} component${plural(count)}`;
    const terms = container.querySelector(".associated-cycle-terms");
    if (!cycle.terms.length) {
      terms.append(element("p", "cycle-empty", "No associated-cycle terms returned."));
      return;
    }

    cycle.terms.forEach((term, index) => {
      const card = element("article", "cycle-term");
      const header = element("div", "cycle-term-header");
      header.append(
        element("span", "cycle-term-number", `Component ${index + 1}`),
        element("span", "cycle-coefficient", `coefficient ${term.coefficient}`),
      );
      card.append(header);

      if (term.underlyingPartition.length) {
        card.append(
          element(
            "span",
            "cycle-partition",
            `underlying partition (${term.underlyingPartition.join(", ")})`,
          ),
        );
      }

      const diagram = element("div", "marked-diagram");
      diagram.setAttribute("role", "img");
      diagram.setAttribute(
        "aria-label",
        term.markedRows.length
          ? `Marked Young diagram with rows ${term.markedRows.join("; ")}.`
          : "Trivial marked Young diagram.",
      );
      if (!term.markedRows.length) {
        diagram.append(element("span", "marked-diagram-empty", "(trivial)"));
      }
      term.markedRows.forEach((row, rowIndex) => {
        const rowElement = element("div", "marked-diagram-row");
        [...row].forEach((marker, columnIndex) => {
          const cell = element(
            "span",
            `marked-diagram-cell cycle-marker-${markedDiagramMarkerClass(marker)}`,
            marker,
          );
          cell.setAttribute("aria-hidden", "true");
          cell.title = `row ${rowIndex + 1}, column ${columnIndex + 1}: ${markedDiagramMarkerLabel(marker)}`;
          rowElement.append(cell);
        });
        diagram.append(rowElement);
      });
      card.append(diagram);
      terms.append(card);
    });
  }

  function markedDiagramMarkerClass(marker) {
    const classes = { "+": "plus", "-": "minus", "*": "star", "=": "equal" };
    return classes[marker] || "other";
  }

  function markedDiagramMarkerLabel(marker) {
    const labels = {
      "+": "trivial local system on a plus sign",
      "-": "trivial local system on a minus sign",
      "*": "nontrivial local system on a plus sign",
      "=": "nontrivial local system on a minus sign",
    };
    return labels[marker] || marker;
  }

  function renderTableau(container, columns, name) {
    container.setAttribute("role", "img");
    container.setAttribute("aria-label", diagramAriaLabel(name, columns));
    if (!columns.length) {
      container.append(element("span", "tableau-empty", "empty diagram"));
      return;
    }

    columns.forEach((column, columnIndex) => {
      const columnElement = element("div", "tableau-column");
      [...column].forEach((marker, rowIndex) => {
        const visibleMarker = marker === "*" ? "•" : marker;
        const cell = element("span", `tableau-cell marker-${markerClass(marker)}`, visibleMarker);
        cell.setAttribute("aria-hidden", "true");
        cell.title = `${name}, column ${columnIndex + 1}, row ${rowIndex + 1}: ${
          marker === "*" ? "* (dot)" : marker
        }`;
        columnElement.append(cell);
      });
      container.append(columnElement);
    });
  }

  function diagramAriaLabel(name, columns) {
    if (!columns.length) return `${name} is an empty diagram.`;
    const spoken = columns.map((column) =>
      [...column].map((marker) => (marker === "*" ? "dot" : marker)).join(" "),
    );
    return `${name} diagram. Columns, left to right: ${spoken.join("; ")}.`;
  }

  function markerClass(marker) {
    const classes = { "*": "star", c: "c", d: "d", s: "s", r: "r" };
    return classes[marker] || "other";
  }

  function renderPath(path, index, k, groupIndex, groupKind) {
    const details = element("details", "path-card");
    details.id = `path-${k}-${groupIndex + 1}-${path.number ?? index + 1}`;

    const summary = document.createElement("summary");
    summary.append(
      element("span", "path-index", path.subsetLabel || String(path.number ?? index + 1)),
    );
    const summaryMain = element("span", "path-summary-main");
    if (path.targetLabel) {
      summaryMain.append(element("code", "path-target", path.targetLabel));
    }
    summaryMain.append(
      element("span", "path-name", path.name),
      element("span", "path-preview", pathHistoryPreview(path.realizations)),
    );
    summary.append(summaryMain);
    const historyCount = path.realizations.length;
    if (groupKind === "so" && path.outerEpsilon) {
      summary.append(element("span", "history-count", `outer ε · ${path.outerEpsilon}`));
    }
    const cycleRelation = path.realizations.find(
      (realization) => realization.pbpCycleRelation,
    )?.pbpCycleRelation;
    if (groupKind === "so" && cycleRelation) {
      const relationClass = cycleRelation === "tensor_1_1" ? "twisted" : "same";
      const relationText = cycleRelation === "tensor_1_1"
        ? "marked cycle · ⊗ (1,1)"
        : "marked cycle · as displayed";
      summary.append(
        element(
          "span",
          `path-cycle-status path-cycle-status-${relationClass}`,
          relationText,
        ),
      );
    }
    summary.append(
      element("span", "history-count", `${historyCount || 1} twist histor${historyCount === 1 ? "y" : "ies"}`),
      element("span", "path-chevron"),
    );
    details.append(summary);

    const content = element("div", "path-content");
    path.realizations.forEach((realization, realizationIndex) => {
      const realizationBlock = element("div", "concrete-realization");
      if (path.realizations.length > 1) {
        realizationBlock.append(
          element("span", "realization-label", `Concrete realization ${realizationIndex + 1}`),
        );
      }
      const histories = element("div", "twist-histories");
      histories.append(element("span", "content-label", "Concrete twist history"));
      const historyList = element("div", "history-list");
      const formatted = realization.history.length
        ? realization.historyText
        : "direct Mp(0) lift";
      historyList.append(element("code", "history", formatted));
      histories.append(historyList);
      realizationBlock.append(histories);

      if (groupKind === "so" && realization.exactAssociatedCycle) {
        realizationBlock.append(renderPathAssociatedCycle(realization));
      }

      realizationBlock.append(element("span", "content-label", "Theta-lift chain"));
      const timeline = element("ol", "chain-timeline");
      if (realization.steps.length) {
        realization.steps.forEach((step) => timeline.append(element("li", "chain-step", step)));
      } else {
        timeline.append(element("li", "chain-step", "No chain steps were supplied."));
      }
      realizationBlock.append(timeline);
      content.append(realizationBlock);
    });
    details.append(content);
    return details;
  }

  function renderPathAssociatedCycle(realization) {
    const panel = element("section", "path-associated-cycle");
    panel.setAttribute("aria-label", "Exact O-path associated cycle");

    const heading = element("header", "path-cycle-heading");
    const title = element("div");
    title.append(
      element("span", "cycle-kicker", "Exact O-path cycle"),
      element("h6", "", "Marked diagrams for this realization"),
    );
    const isTwisted = realization.pbpCycleRelation === "tensor_1_1";
    heading.append(
      title,
      element(
        "span",
        `path-cycle-relation path-cycle-relation-${isTwisted ? "twisted" : "same"}`,
        isTwisted ? "reference ⊗ (1,1)" : "same as reference",
      ),
    );
    panel.append(heading);
    panel.append(
      element(
        "p",
        "path-cycle-explanation",
        isTwisted
          ? "This exact path reaches the painted-bipartition cycle above twisted by ⊗ (1,1)."
          : "This exact path reaches the painted-bipartition cycle displayed above.",
      ),
    );

    const cycle = element("div", "path-cycle-render");
    cycle.append(
      element("span", "cycle-count"),
      element("div", "associated-cycle-terms"),
    );
    panel.append(cycle);
    renderAssociatedCycle(cycle, realization.exactAssociatedCycle);
    return panel;
  }

  function updateToggleLabel(button, details) {
    if (!details.length) {
      button.hidden = true;
      return;
    }
    button.textContent = details.every((detail) => detail.open) ? "Close all" : "Open all";
  }

  function pathHistoryPreview(realizations) {
    const histories = realizations
      .map((realization) => realization.historyText)
      .filter(Boolean);
    if (!histories.length) return "Direct Mp(0) lift";
    return histories.join(" · ");
  }

  function uniquePathCount(groups) {
    const ids = new Set();
    let anonymous = 0;
    groups.forEach((group) => {
      group.paths.forEach((path) => {
        if (path.id === undefined || path.id === null) anonymous += 1;
        else ids.add(String(path.id));
      });
    });
    return ids.size + anonymous;
  }

  function toNumberArray(value) {
    return Array.isArray(value) ? value.map(Number).filter(Number.isFinite) : [];
  }

  function stringArray(value) {
    if (!Array.isArray(value)) return [];
    return value.map(String);
  }

  function numberOr(value, fallback) {
    const number = Number(value);
    return Number.isFinite(number) ? number : fallback;
  }

  function plural(number) {
    return number === 1 ? "" : "s";
  }

  function strongText(value) {
    return element("strong", "", String(value));
  }

  function element(tag, className = "", text = "") {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== "") node.textContent = text;
    return node;
  }

  function reducedMotion() {
    return window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  }
})();
