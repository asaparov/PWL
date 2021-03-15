document.addEventListener('DOMContentLoaded', function(event) {
	const socket = io();

	var console = document.getElementById('console');
	var input = document.getElementById('input');

	input.addEventListener('keyup', function(event) {
		if (event.keyCode === 13) {
			event.preventDefault();
			socket.emit('input', input.value);
			input.value = '';
			console.innerHTML += '\n';
		}
	});

	socket.on('console', (data) => {
		console.innerHTML += data
			.replace(/&/g, "&amp;")
			.replace(/</g, "&lt;")
			.replace(/>/g, "&gt;")
			.replace(/"/g, "&quot;")
			.replace(/'/g, "&#039;");
	});
});
