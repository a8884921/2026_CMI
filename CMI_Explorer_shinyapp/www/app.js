/* CMI Explorer: page navigation (Analysis with its three workspaces, Help, About; history, page
   titles, deep links, mobile menu), accessibility helpers (skip link, main landmark, tab semantics)
   and the small messages the server sends (enable a button, move focus, scroll, warn before leaving
   during an analysis). */
(function () {
  'use strict';
  var INNER = ['single_cell', 'spatial', 'bulk'];          // workspaces inside the Analysis page
  var busy = {}, fromPop = false, pendingTarget = null, initial = true;

  function navLink(v) {
    return document.querySelector('.navbar a[data-value="' + v + '"][data-bs-toggle], #ws a[data-value="' + v + '"][data-bs-toggle]');
  }
  function show(link) {
    if (window.bootstrap && bootstrap.Tab) bootstrap.Tab.getOrCreateInstance(link).show();
    else if (window.jQuery) jQuery(link).tab('show');
  }
  function paneOf(link) { var h = link.getAttribute('href') || link.getAttribute('data-bs-target') || ''; return h.charAt(0) === '#' ? document.getElementById(h.slice(1)) : null; }
  function activeOuter() { var a = document.querySelector('.navbar a.active[data-value][data-bs-toggle]'); return a ? a.getAttribute('data-value') : null; }
  function activeInner() { var a = document.querySelector('#ws a.active[data-value]'); return a ? a.getAttribute('data-value') : INNER[0]; }
  function focusEl(el) { if (!el.hasAttribute('tabindex') && !/^(A|BUTTON|INPUT|SELECT|TEXTAREA)$/.test(el.tagName)) el.setAttribute('tabindex', '-1'); el.focus({ preventScroll: true }); }

  // Show a page: analysis, help, about, or one of the workspaces (which opens Analysis).
  function showPage(v) {
    var outer = INNER.indexOf(v) >= 0 ? 'analysis' : v, o = navLink(outer);
    if (!o) return false;
    var moved = false;
    if (outer !== v) { var i = navLink(v); if (i && !i.classList.contains('active')) { show(i); moved = true; } }
    if (!o.classList.contains('active')) { show(o); return true; }
    if (!moved) settle();
    return true;
  }

  // After a page is shown: current-page state, title, history, mobile menu, scroll and focus.
  function settle() {
    var outer = activeOuter(), page = outer === 'analysis' ? activeInner() : outer, o = navLink(outer), p = navLink(page) || o;
    document.querySelectorAll('.navbar [aria-current="page"]').forEach(function (el) { el.removeAttribute('aria-current'); });
    if (o) o.setAttribute('aria-current', 'page');
    document.title = (p ? p.textContent.trim() + ' · ' : '') + 'CMI Explorer';
    if (!fromPop && page && window.history && history.pushState && location.hash !== '#' + page) {
      if (initial) history.replaceState({ tab: page }, '', '#' + page); else history.pushState({ tab: page }, '', '#' + page);
    }
    fromPop = false; initial = false;
    var open = document.querySelector('.navbar-collapse.show');
    if (open && window.bootstrap && bootstrap.Collapse) bootstrap.Collapse.getOrCreateInstance(open).hide();
    var target = pendingTarget; pendingTarget = null;
    setTimeout(function () {
      if (target && target.el) {
        if (target.el.tagName === 'DETAILS') target.el.open = true;
        target.el.scrollIntoView({ block: 'start' }); if (target.focus) focusEl(target.el); return;
      }
      window.scrollTo(0, 0);
      if (o && o.getAttribute('data-user') === '1') {
        var pane = paneOf(o), h = pane ? pane.querySelector('h1') : null;
        if (h) { h.setAttribute('tabindex', '-1'); h.focus({ preventScroll: true }); }
        o.removeAttribute('data-user');
      }
    }, 30);
  }

  // #help, #single_cell ... open a page; #h-spatial (any id inside a page) opens it and scrolls there.
  function pageOf(el) {
    var inner = null, outer = null, pane = el.closest('.tab-pane');
    while (pane) {
      var link = document.querySelector('.navbar a[data-bs-toggle][href="#' + pane.id + '"], #ws a[data-bs-toggle][href="#' + pane.id + '"]');
      if (link) { var v = link.getAttribute('data-value'); if (INNER.indexOf(v) >= 0) inner = inner || v; else outer = outer || v; }
      pane = pane.parentElement ? pane.parentElement.closest('.tab-pane') : null;
    }
    return inner || outer;
  }
  function fromHash() {
    var h = decodeURIComponent((location.hash || '').replace(/^#/, ''));
    if (!h) return showPage('analysis');
    if (navLink(h)) return showPage(h);
    var el = document.getElementById(h), v = el ? pageOf(el) : null;
    if (!v) return false;
    pendingTarget = { el: el, focus: false };
    return showPage(v);
  }

  // Result tabs (bslib navset_tab inside the workspaces): tab roles, selected state and arrow keys.
  function ariaTabs(root) {
    (root || document).querySelectorAll('.result-tabs .nav-tabs').forEach(function (list) {
      list.setAttribute('role', 'tablist');
      list.querySelectorAll('a.nav-link').forEach(function (a, i) {
        var pane = paneOf(a);
        if (!a.id) a.id = (pane ? pane.id : 'restab' + i) + '-tab';
        a.setAttribute('role', 'tab');
        var sel = a.classList.contains('active');
        a.setAttribute('aria-selected', sel ? 'true' : 'false');
        a.setAttribute('tabindex', sel ? '0' : '-1');
        if (pane) { a.setAttribute('aria-controls', pane.id); pane.setAttribute('role', 'tabpanel'); pane.setAttribute('aria-labelledby', a.id); pane.setAttribute('tabindex', '0'); }
        if (a.parentElement) a.parentElement.setAttribute('role', 'presentation');
      });
      if (list.getAttribute('data-keys')) return;
      list.setAttribute('data-keys', '1');
      list.addEventListener('keydown', function (e) {
        var tabs = Array.prototype.slice.call(list.querySelectorAll('a.nav-link'));
        var i = tabs.indexOf(document.activeElement); if (i < 0) return;
        var j = e.key === 'ArrowRight' ? (i + 1) % tabs.length : e.key === 'ArrowLeft' ? (i - 1 + tabs.length) % tabs.length : e.key === 'Home' ? 0 : e.key === 'End' ? tabs.length - 1 : -1;
        if (j < 0) return;
        e.preventDefault(); tabs[j].focus(); show(tabs[j]);
      });
    });
  }

  document.addEventListener('shown.bs.tab', function (e) {
    var t = e.target;
    if (!t.closest) return;
    if (t.closest('.navbar')) { settle(); return; }
    if (t.closest('#ws')) { if (activeOuter() === 'analysis') settle(); return; }
    if (t.closest('.result-tabs')) ariaTabs(t.closest('.result-tabs').parentElement);
  });
  window.addEventListener('popstate', function () { fromPop = true; if (!fromHash()) fromPop = false; });

  document.addEventListener('DOMContentLoaded', function () {
    // one main landmark: the container of the pages
    var tc = document.querySelector('body > .container-fluid > .tab-content') || document.querySelector('.tab-content');
    if (tc) { tc.setAttribute('role', 'main'); tc.id = tc.id || 'main-content'; }
    var skip = document.createElement('a');
    skip.className = 'skip-link'; skip.href = '#main-content'; skip.textContent = 'Skip to main content';
    skip.addEventListener('click', function (e) {
      e.preventDefault();
      var o = navLink(activeOuter() || 'analysis'), p = o ? paneOf(o) : null, h = p ? p.querySelector('h1') : null;
      if (h) { h.setAttribute('tabindex', '-1'); h.focus(); }
    });
    document.body.insertBefore(skip, document.body.firstChild);
    document.querySelectorAll('.navbar a[data-value][data-bs-toggle]').forEach(function (a) {
      a.addEventListener('click', function () {
        a.setAttribute('data-user', '1');
        var open = document.querySelector('.navbar-collapse.show');            // also when the page is already shown
        if (open && window.bootstrap && bootstrap.Collapse) bootstrap.Collapse.getOrCreateInstance(open).hide();
      });
    });
    if (!fromHash()) settle();

    document.body.addEventListener('click', function (e) {
      var go = e.target.closest ? e.target.closest('[data-go]') : null;
      if (!go) return;
      e.preventDefault();
      var page = go.getAttribute('data-go'), sub = go.getAttribute('data-subtab'), ex = go.getAttribute('data-example');
      if (sub) { var el = document.getElementById(sub); if (el) pendingTarget = { el: el, focus: false }; }
      showPage(page);
      if (ex) setTimeout(function () { var b = document.getElementById(ex); if (b && !b.disabled) b.click(); }, 200);
    });

    // Result panels are rendered by the server; give their tabs roles as they appear.
    if (window.MutationObserver) new MutationObserver(function (muts) {
      for (var i = 0; i < muts.length; i++) for (var j = 0; j < muts[i].addedNodes.length; j++) {
        var n = muts[i].addedNodes[j];
        if (n.nodeType === 1 && (n.matches('.result-tabs') || n.querySelector('.result-tabs'))) { ariaTabs(n.parentElement || n); return; }
      }
    }).observe(document.body, { childList: true, subtree: true });
  });

  // A new file is not checked yet: lock "Run analysis" until the server's check of this upload arrives.
  document.addEventListener('change', function (e) {
    var t = e.target;
    if (!t || t.type !== 'file' || !t.closest || !t.closest('.ws-side')) return;
    var m = (t.id || '').match(/^(.*)-(file|clinical)$/); if (!m) return;
    if (!t.files || !t.files.length) return;   // dialog cancelled: the previous upload and its check stand
    var run = document.getElementById(m[1] + '-run');
    if (run) { run.disabled = true; run.setAttribute('aria-disabled', 'true'); }
    var f = t.files[0], maxMb = parseFloat(t.getAttribute('data-max-mb') || '');
    if (maxMb > 0 && f.size > maxMb * 1048576 && window.Shiny && Shiny.setInputValue) {
      // Above the server's upload limit: do not send it; the server records the refusal and explains it.
      e.stopPropagation();
      var group = t.closest('.input-group'), txt = group && group.querySelector('input[type=text]');
      if (txt) txt.value = f.name;
      var prog = document.getElementById(t.id + '_progress'); if (prog) prog.style.visibility = 'hidden';   // not "Upload complete" of the earlier file
      t.value = '';
      Shiny.setInputValue(t.id + '_refused', { name: f.name, mb: f.size / 1048576 }, { priority: 'event' });
      return;
    }
    var box = document.getElementById(m[1] + (m[2] === 'file' ? '-file_check' : '-clinical_check'));
    if (box) box.innerHTML = '<div class="file-check checking"><b>Uploading and checking the file…</b><span>The analysis can start when the check has passed.</span></div>';
  }, true);

  // Scrolling table bodies of DataTables can be reached with the keyboard (WCAG 2.1.1).
  if (window.jQuery) jQuery(document).on('init.dt', function (e, settings) {
    var w = settings && settings.nTableWrapper; if (!w) return;
    w.querySelectorAll('.dataTables_scrollBody').forEach(function (b) {
      b.setAttribute('tabindex', '0'); b.setAttribute('role', 'region'); b.setAttribute('aria-label', 'Table, scrolls horizontally');
    });
  });

  window.addEventListener('beforeunload', function (e) {
    for (var k in busy) if (busy[k]) { e.preventDefault(); e.returnValue = ''; return ''; }
  });

  function onShiny() {
    if (!window.Shiny || !Shiny.addCustomMessageHandler || window.__cmiHandlers) return;
    window.__cmiHandlers = true;
    Shiny.addCustomMessageHandler('cmi-enable', function (m) {
      var b = document.getElementById(m.id); if (!b) return;
      b.disabled = !m.enabled; if (m.enabled) b.removeAttribute('aria-disabled'); else b.setAttribute('aria-disabled', 'true');
    });
    Shiny.addCustomMessageHandler('cmi-focus', function (m) {
      var tries = 0;
      (function attempt() {
        var el = document.getElementById(m.id);
        if (el) { focusEl(el); el.scrollIntoView({ block: 'start' }); return; }
        if (++tries < 40) setTimeout(attempt, 100);
      })();
    });
    Shiny.addCustomMessageHandler('cmi-scroll', function (m) {
      var el = document.getElementById(m.id); if (!el) return;
      var r = el.getBoundingClientRect();
      if (r.top < 0 || r.top > window.innerHeight * 0.6) el.scrollIntoView({ block: 'start' });
    });
    Shiny.addCustomMessageHandler('cmi-busy', function (m) { busy[m.id] = !!m.busy; });
    setInterval(function () { if (Shiny.setInputValue) Shiny.setInputValue('keepalive', Date.now(), { priority: 'event' }); }, 240000);
  }
  onShiny();
  document.addEventListener('DOMContentLoaded', onShiny);
  if (window.jQuery) jQuery(document).on('shiny:connected', onShiny);
})();
