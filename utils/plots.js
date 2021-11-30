

/**
 * Plots a data set to the document
 * 
 * @param {array} x Set of labels to show
 * @param {array} y Data set to show
 * @param {string} name Label to attach to the data set
 * @param {string} colour HTML colour to use for the line
 */
function plot (x, y, name, colour = 'black') {

    // Create new 2D canvas to hold the results
    let canvas = document.createElement('canvas');
    canvas.width = window.innerWidth * 0.6;
    canvas.height = window.innerHeight * 0.6;
    canvas.style.margin = '50px auto';
    canvas.style.display = 'block';
    canvas.style.background = 'white';
    canvas.style.border = 'solid black 1px';
    canvas.style.padding = '10px';
    document.body.append(canvas);

    // Plot through Chart.js
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

