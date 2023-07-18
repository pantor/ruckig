import RuckigModule from './ruckig.js';
import Alpine from 'alpinejs';

var Plotly = require('plotly.js/lib/core');


RuckigModule().then(lib => {
    function toWASM(vec) {
        const q = new lib.Vector();
        q.resize(vec.length, 0.0);
        for (let i = 0; i < vec.length; i += 1) q.set(i, vec[i]);
        return q;
    }

    Alpine.data('dropdown', () => ({
        current_position: 0.0,
        current_velocity: 0.0,
        current_acceleration: 0.0,
        target_position: 1.0,
        target_velocity: 0.0,
        target_acceleration: 0.0,
        max_velocity: 1.0,
        max_acceleration: 1.0,
        max_jerk: 2.0,

        updatePlot() {
            const input = new lib.InputParameter(1);
            input.current_position = toWASM([this.current_position]);
            input.current_velocity = toWASM([this.current_velocity]);
            input.current_acceleration = toWASM([this.current_acceleration]);
            input.target_position = toWASM([this.target_position]);
            input.target_velocity = toWASM([this.target_velocity]);
            input.target_acceleration = toWASM([this.target_acceleration]);
            input.max_velocity = toWASM([this.max_velocity]);
            input.max_acceleration = toWASM([this.max_acceleration]);
            input.max_jerk = toWASM([this.max_jerk]);

            const trajectory = new lib.Trajectory(1);

            const ruckig = new lib.Ruckig(1);
            const result = ruckig.calculate(input, trajectory);

            const duration = trajectory.get_duration();

            if (result.value !== 0) {
                console.log('error');
                return;
            }

            let ts = [], ps = [], vs = [], as = [], js = [];
            for (let t = 0; t < duration; t += duration / 100) {
                const state = trajectory.at_time(t);

                ts.push(t);
                ps.push(state.position.get(0));
                vs.push(state.velocity.get(0));
                as.push(state.acceleration.get(0));
                js.push(state.jerk.get(0));
            }
            const data = [
                {x: ts, y: ps, name: 'Position'},
                {x: ts, y: vs, name: 'Velocity'},
                {x: ts, y: as, name: 'Acceleration'},
                {x: ts, y: js, name: 'Jerk'},
            ];
            const layout = {
                title: 'Trajectory',
                    autosize: true,
                    dragmode: 'pan',
                    margin: {
                        l: 50, r: 10,
                        b: 70, t: 50,
                        pad: 10
                    },
                hovermode: "x",
                shapes: [{
                    type: 'line',
                    xref: 'x', yref: 'paper',
                    x0: duration, y0: -10,
                    x1: duration, y1: 10,
                    line: {
                        color: 'black',
                        width: 0.5,
                    }
                }]
            };
            const config = {
                modeBarButtonsToRemove: ['zoom'],
                scrollZoom: true,
            }
            Plotly.newPlot('plot', data, layout, config);
        }
    }));

    Alpine.start();
});
