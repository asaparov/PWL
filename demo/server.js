const express = require('express');
const app = express();
const path = require('path');
const server = require('http').createServer(app);
const io = require("socket.io")(server);
const spawn = require('child_process').spawn;
const port = 3000;

server.listen(port, () => {
	console.log('Server listening at port %d.', port);
});

app.use(express.static(__dirname));

var id = 0;

io.on("connection", socket => {
	var current_id = id;
	id++;
	console.log('New connection with ID %d.', current_id);
	const process = spawn('/usr0/home/asaparov/theory_induction/executive_test_dbg', [], {cwd: '/usr0/home/asaparov/theory_induction'});

	process.stdout.on('data', (data) => {
		socket.emit('console', data.toString());
	});

	process.stderr.on('data', (data) => {
		socket.emit('console', data.toString());
	});

	process.on('close', (code) => {
		socket.emit('console', 'Application crashed!');
	});

	socket.on('input', (data) => {
		process.stdin.write(data + '\n');
	});

	socket.on('disconnect', () => {
		console.log('Connection %d disconnected.', current_id);
		process.kill();
	});
});
