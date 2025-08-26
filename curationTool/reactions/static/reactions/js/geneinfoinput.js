function createGeneInfoInput() {
	const wrapper = document.getElementById('geneInfoWrapper');
	if (!wrapper) return;

	// --- Prevent duplicate init (double textareas & stale handlers) ---
	if (document.getElementById('geneInfoInput')) return;

	const container = document.createElement('div');
	container.id = 'geneInfoContainer';
	container.className = 'gene-info-container';
	container.style.position = 'relative';

	const ta = document.createElement('textarea');
	ta.id = 'geneInfoInput';
	ta.className = 'gene-info-input';
	ta.rows = 1;
	ta.placeholder = 'Type HGNC symbol or Entrez ID. Use AND / OR / ( ) to build GPR';
	Object.assign(ta.style, {
		width: '100%',
		padding: '10px',
		border: '1px solid #ccc',
		borderRadius: '6px',
		fontFamily: 'inherit',
		resize: 'vertical',
	});

	const box = document.createElement('div');
	box.id = 'geneSuggest';
	box.className = 'suggest hidden';
	Object.assign(box.style, {
		position: 'absolute',
		left: 0,
		right: 0,
		background: '#fff',
		border: '1px solid rgba(34,36,38,.15)',
		borderTop: 'none',
		boxShadow: '0 10px 30px rgba(0,0,0,.1)',
		zIndex: 10,
		maxHeight: '260px',
		overflowY: 'auto',
	});

	container.appendChild(ta);
	container.appendChild(box);
	wrapper.appendChild(container);

	// --- ensure dropdown renders below the textarea ---
	function positionBox() {
		// place just under the textarea inside the relatively positioned container
		box.style.top = ta.offsetHeight + 'px';
	}
	positionBox();
	const ro = new ResizeObserver(positionBox);
	ro.observe(ta);
	window.addEventListener('resize', positionBox);

	let items = [];
	let active = -1;
	let debounceId = null;
	let lastTokenInfo = null;
	let navigating = false; // true after ArrowUp/Down
	let picking = false; // true while clicking inside dropdown

	function getTokenBounds(value, caret) {
		const isWord = (ch) => /[A-Za-z0-9]/.test(ch);
		let s = caret,
			e = caret;
		while (s > 0 && isWord(value[s - 1])) s--;
		while (e < value.length && isWord(value[e])) e++;
		const token = value.slice(s, e);
		return { start: s, end: e, token };
	}

	async function fetchSuggest(q) {
		try {
			const res = await fetch(`/api/gene/suggest?q=${encodeURIComponent(q)}`);
			const data = await res.json();
			return data.items || [];
		} catch {
			return [];
		}
	}

	function renderSuggest() {
		if (!items.length) {
			box.classList.add('hidden');
			box.innerHTML = '';
			return;
		}
		box.innerHTML = items
			.map(
				(it, i) => `
      <div class="s-item ${i === active ? 'active' : ''}" data-i="${i}"
           style="display:flex;gap:.5rem;padding:.5rem .75rem;cursor:pointer;">
        <div class="s-sym" style="font-weight:700;">${it.symbol}</div>
        <div class="s-name" style="opacity:.8;flex:1;">${it.name || ''}</div>
        <div class="s-badge ${it.present ? 'ok' : 'new'}"
             style="font-size:.75rem;opacity:.75;white-space:nowrap;">
          ${it.present ? 'in VMH' : 'new'}
        </div>
      </div>
    `
			)
			.join('');
		box.classList.remove('hidden');
		positionBox();
	}

	// Accept only when explicitly chosen
	function acceptItem(it) {
		const symbol = it.symbol;
		const v = ta.value;
		const { start, end } = lastTokenInfo || getTokenBounds(v, ta.selectionStart);
		ta.value = v.slice(0, start) + symbol + v.slice(end);
		const caret = start + symbol.length;
		ta.setSelectionRange(caret, caret);
		items = [];
		active = -1;
		navigating = false;
		renderSuggest();
	}

	// Only accept on explicit CLICK inside the dropdown
	box.addEventListener('mousedown', () => {
		picking = true;
	}); // avoid blur closing before click
	box.addEventListener('click', (e) => {
		const row = e.target.closest('.s-item');
		if (!row) return;
		const idx = +row.dataset.i;
		acceptItem(items[idx]);
		picking = false;
		ta.focus();
	});

	// Close suggestions on ANY outside click before it becomes a click
	document.addEventListener(
		'mousedown',
		(e) => {
			if (box.classList.contains('hidden')) return;
			if (e.target === ta || box.contains(e.target)) return;
			items = [];
			active = -1;
			navigating = false;
			renderSuggest();
		},
		true
	); // capture phase = before target handlers

	ta.addEventListener('input', () => {
		const v = ta.value;
		const caret = ta.selectionStart;
		const { start, end, token } = getTokenBounds(v, caret);
		lastTokenInfo = { start, end, token, mode: /^\d+$/.test(token) ? 'Entrez ID' : 'HGNC Symbol' };
		navigating = false;

		const T = token.toUpperCase();
		if (!token || T === 'AND' || T === 'OR') {
			items = [];
			renderSuggest();
			return;
		}
		if (lastTokenInfo.mode === 'HGNC Symbol' && token.length < 2) {
			items = [];
			renderSuggest();
			return;
		}

		clearTimeout(debounceId);
		debounceId = setTimeout(async () => {
			items = await fetchSuggest(token);
			active = items.length ? 0 : -1;
			renderSuggest();
		}, 180);
	});

	ta.addEventListener('keydown', (e) => {
		if (box.classList.contains('hidden')) return;

		if (e.key === 'ArrowDown') {
			e.preventDefault();
			navigating = true;
			active = items.length ? (active + 1) % items.length : -1;
			renderSuggest();
		} else if (e.key === 'ArrowUp') {
			e.preventDefault();
			navigating = true;
			active = items.length ? (active - 1 + items.length) % items.length : -1;
			renderSuggest();
		} else if (e.key === 'Enter' || e.key === 'Tab') {
			// ONLY accept if user navigated suggestions
			if (items.length && navigating) {
				e.preventDefault();
				acceptItem(items[active >= 0 ? active : 0]);
			}
			// otherwise: do nothing special
		} else if (e.key === 'Escape') {
			items = [];
			active = -1;
			navigating = false;
			renderSuggest();
		}
	});

	ta.addEventListener('blur', () => {
		setTimeout(() => {
			if (!picking) {
				items = [];
				active = -1;
				navigating = false;
				renderSuggest();
			}
			picking = false;
		}, 0);
	});

	// Helper the submit code can call before reading the value
	window.closeGeneSuggestions = function () {
		items = [];
		active = -1;
		navigating = false;
		renderSuggest();
	};
}
