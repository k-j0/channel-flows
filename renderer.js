
/**
 * Simple 2D renderer to output a dynamic texture to the screen
 */
class Renderer {

    width; // width of the texture
    height; // height of the texture
    data; // texture data being rendered


    /**
     * Default constructor - initializes the renderer
     * 
     * @param {int} texWidth The width of the texture to display
     * @param {int} texHeight The height of the texture to display
     */
    constructor(texWidth, texHeight) {

        this.width = texWidth;
        this.height = texHeight;

        // Initialise Three.js renderer
        const width = window.innerWidth * 0.9;
        const height = width * 0.2;
        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(width, height);
        renderer.domElement.style.margin = '20px auto';
        document.body.appendChild(renderer.domElement);
        const camera = new THREE.OrthographicCamera(-0.5, 0.5, 0.5, -0.5, 0.1, 10); // x=-0.5..0.5, y=-0.5..0.5 view

        // Create data texture to hold results
        this.data = new Uint8Array(3 * texWidth * texHeight);
        const texture = new THREE.DataTexture(this.data, texWidth, texHeight, THREE.RGBFormat);

        // Show texture on 0..1 plane mesh
        const scene = new THREE.Scene();
        const plane = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), new THREE.MeshBasicMaterial({
            map: texture // bind data texture to plane mesh
        }));
        plane.position.z = -1;
        scene.add(plane);

        // Render function
        this.refresh = () => {
            texture.needsUpdate = true;
            renderer.render(scene, camera);
        };
    }

    /**
     * Refreshes the render by updating the texture
     */
    render(updateCallback) {
        // Gather values to show
        let tmp = [];
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {
                let val = updateCallback(x, y) * 255;
                tmp.push(val);
                if(val < min) min = val;
                if(val > max) max = val;
            }
        }
        // Update texture data buffer with remapped values, mapped to a nice blue colour gradient
        let colours = [
            new THREE.Color('#001D3A'),
            new THREE.Color('#5A8AA6'),
            new THREE.Color('#E3FFFF')
        ];
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {
                let idx = (x + y * this.width) * 3;
                let val = tmp.shift();
                val = (val - min) / (max - min); // remap to 0..1 no matter what
                let col = colours[colours.length - 1];
                // apply colour map
                if (val < 1.0) {
                    val *= colours.length - 1;
                    let colourIdx = Math.floor(val);
                    col = (new THREE.Color()).lerpColors(colours[colourIdx], colours[colourIdx+1], val - colourIdx);
                }
                this.data[idx] = col.r * 255;
                this.data[idx+1] = col.g * 255;
                this.data[idx+2] = col.b * 255;
            }
        }
        this.refresh();
    }


}; // class Renderer
