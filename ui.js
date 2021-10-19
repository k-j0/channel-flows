
/**
 * Shows a button on the document with a bound callback
 * 
 * @param {string} name Text to display on the button itself
 * @param {function} callback Function to call upon clicks being registered
 */
function button(name, callback) {
    let button = document.createElement('button');
    button.innerHTML = name;
    button.style.margin = '10px auto';
    button.style.display = 'block';
    button.addEventListener('click', callback);
    document.body.append(button);
}
