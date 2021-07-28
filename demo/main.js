function hsv_to_rgb(h, s, v) {
	h_i = parseInt(h*6);
	f = h*6 - h_i;
	p = v * (1 - s);
	q = v * (1 - f*s);
	t = v * (1 - (1 - f) * s);
	switch (h_i) {
	case 0: rgb = [v, t, p]; break;
	case 1: rgb = [q, v, p]; break;
	case 2: rgb = [p, v, t]; break;
	case 3: rgb = [p, q, v]; break;
	case 4: rgb = [t, p, v]; break;
	case 5: rgb = [v, p, q]; break;
	}
	return [parseInt(rgb[0]*256), parseInt(rgb[1]*256), parseInt(rgb[2]*256)];
}

function highlight(x) {
	e = document.getElementsByClassName(x.className);
	for (i = 0; i < e.length; i++) {
		el = e[i];
		el.style.textShadow = '0 0 25px #000, 0 0 5px #2272d3';
	}
}

function unhighlight(x) {
	e = document.getElementsByClassName(x.className);
	for (i = 0; i < e.length; i++) {
		el = e[i];
		el.style.textShadow = 'none';
	}
}

const phi_conjugate = 0.618033988749895;
var h = 0.0;
var next_color = 0;

document.addEventListener('DOMContentLoaded', function(event) {
	const socket = io();
	const CSI   = "\u001b[";
	const BLUE  = CSI + "34m";
	const BOLD  = CSI + "1m";
	const RESET = CSI + "0m";
	const SUBSCRIPTS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];

	var console = document.getElementById('console');
	var input = document.getElementById('input');

	input.addEventListener('keyup', function(event) {
		if (event.keyCode === 13) {
			event.preventDefault();
			socket.emit('input', input.value);
			console.innerHTML += input.value + '\n';
			input.value = '';
		}
	});

	socket.on('console', (data) => {
		var new_data = [];
		var style_stack = [];
		var sub = false;
		for (var i = 0; i < data.length; i++) {
			if (data[i] == '&') {
				new_data.push("&amp;");
			} else if (data[i] == '<') {
				new_data.push("&lt;");
			} else if (data[i] == '>') {
				new_data.push("&gt;");
			} else if (data[i] == '"') {
				new_data.push("&quot;");
			} else if (data[i] == '\'') {
				new_data.push("&#039;");
			} else if (data.startsWith(BOLD, i)) {
				new_data.push("<b>");
				style_stack.push("b");
				i += BOLD.length - 1;
			} else if (data.startsWith(BLUE, i)) {
				new_data.push("<span class=\"blue\">");
				style_stack.push("span");
				i += BLUE.length - 1;
			} else if (data.startsWith(RESET, i)) {
				while (style_stack.length != 0)
					new_data.push("</" + style_stack.pop() + ">");
				i += RESET.length - 1;
			} else if ((data[i] == 'x' || data[i] == 'c') && i + 1 < data.length && SUBSCRIPTS.includes(data[i + 1])) {
				j = i + 1;
				var index = SUBSCRIPTS.indexOf(data[j]);
				var digits = [];
				do {
					digits.push(index.toString());
					index = SUBSCRIPTS.indexOf(data[++j]);
				} while (index != -1);
				const num_str = digits.join('');
				const num = parseInt(num_str);
				if (next_color <= num) {
					var new_style = ["<style>"];
					while (next_color <= num) {
						[r,g,b] = hsv_to_rgb(h,0.85,0.5);
						new_style.push(".var-" + next_color + "{color:rgb(" + r + "," + g + "," + b + ");border-radius:2px;}");
						h = (h + phi_conjugate) % 1;
						next_color++;
					}
					new_style.push("</style>");
					console.innerHTML += new_style.join("");
				}
				new_data.push("<span class=\"var-" + num_str + "\" onmouseover=\"highlight(this)\" onmouseout=\"unhighlight(this)\">" + data[i] + "<sub>" + num_str + "</sub></span>");
				i = j - 1;
			} else if ((index = SUBSCRIPTS.indexOf(data[i])) != -1) {
				if (!sub) {
					new_data.push("<sub>");
					sub = true;
				}
				new_data.push('' + index);
			} else {
				if (sub) {
					new_data.push("</sub>");
					sub = false;
				}
				new_data.push(data[i]);
			}
		}
		console.innerHTML += new_data.join("");
	});
});
