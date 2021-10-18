

function plot (x, y, name, colour = 'black') {

    let canvas = document.createElement('canvas');
    canvas.width = window.innerWidth * 0.6;
    canvas.height = window.innerHeight * 0.6;
    canvas.style.margin = '50px auto';
    canvas.style.display = 'block';
    canvas.style.background = 'white';
    canvas.style.border = 'solid black 1px';
    canvas.style.padding = '10px';
    document.body.append(canvas);

    let ctx = canvas.getContext('2d');
    new Chart(ctx, {
        type: 'line',
        data: {
            labels: x,
            datasets: [
                {
                    label: name,
                    data: y,
                    borderColor: colour
                }
            ]
        }
    });

}

