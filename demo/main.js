document.addEventListener('DOMContentLoaded', function(event) {
	const socket = io();
	const CSI   = "\u001b[";
	const BLUE  = CSI + "34m";
	const BOLD  = CSI + "1m";
	const RESET = CSI + "0m";

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
			} else {
				new_data.push(data[i]);
			}
		}
		console.innerHTML += new_data.join("");
	});
});
