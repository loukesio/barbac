(function () {
  'use strict';
  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('studio-busy', function (busy) {
      ['run', 'extract', 'build_report', 'reset', 'load_demo'].forEach(function (id) {
        const button = document.getElementById(id);
        if (button) button.disabled = busy;
      });
      document.querySelectorAll('input[type="file"]').forEach(function (el) { el.disabled = busy; });
      document.body.setAttribute('aria-busy', String(busy));
    });
    Shiny.addCustomMessageHandler('studio-clear-files', function (message) {
      document.querySelectorAll('input[type="file"]').forEach(function (el) {
        el.value = '';
        const text = el.closest('.shiny-input-container').querySelector('input[type="text"]');
        if (text) text.value = '';
        Shiny.setInputValue(el.id, null, { priority: 'event' });
      });
    });
    Shiny.addCustomMessageHandler('studio-download-state', function (state) {
      const ids = ['download_all', 'download_centroids', 'download_memberships', 'download_series', 'download_plot'];
      ids.concat(['download_input']).forEach(function (id) {
        const el = document.getElementById(id);
        if (!el) return;
        const ready = id === 'download_input' ? state.input : id === 'download_plot' ? state.plot : state.results;
        el.classList.toggle('disabled', !ready);
        el.setAttribute('aria-disabled', String(!ready));
        el.tabIndex = ready ? 0 : -1;
      });
    });
  });
  $(document).on('change', '#page input', function () {
    window.scrollTo({ top: 0, behavior: 'instant' });
  });
})();
