
/**
 * D2Q9 lattice boltzmann method simulation on rectangular lattice
 */
class D2Q9 {

    /// Simulation constants (always the same across any D2Q9 sim)
    static W = [ // weights
        4.0 / 9.0,  // i = 0
        1.0 / 9.0,  // i = 1
        1.0 / 9.0,  // i = 2
        1.0 / 9.0,  // i = 3
        1.0 / 9.0,  // i = 4
        1.0 / 36.0, // i = 5
        1.0 / 36.0, // i = 6
        1.0 / 36.0, // i = 7
        1.0 / 36.0  // i = 8
    ];
    static C = [ // velocities
        new THREE.Vector2(0, 0),    // i = 0
        new THREE.Vector2(1, 0),    // i = 1
        new THREE.Vector2(0, -1),   // i = 2
        new THREE.Vector2(-1, 0),   // i = 3
        new THREE.Vector2(0, 1),    // i = 4
        new THREE.Vector2(1, -1),   // i = 5
        new THREE.Vector2(-1, -1),  // i = 6
        new THREE.Vector2(-1, 1),   // i = 7
        new THREE.Vector2(1, 1)     // i = 8
    ];
    

    /**
     * Computes the equilibrium value Ni^eq for a site
     * 
     * @param {number} i Direction for which to get Ni^eq
     * @param {number} density The local density ρ
     * @param {THREE.Vector2} velocity The velocity u
     */
    static NiEq (i, density, velocity) {
        let uDotU = velocity.lengthSq(); // u • u
        let ciDotU = D2Q9.C[i].dot(velocity); // c_i • u
        return D2Q9.W[i] * density * (1.0 + 3.0 * ciDotU + 4.5 * ciDotU * ciDotU - 1.5 * uDotU);
    }




    /**
     * Parameters given to initialize the simulation
     */
    static Parameters = class {
        width = 20; // int, number of sites across the lattice horizontally
        height = 20; // int, number of sites across the lattice vertically
        kinematicViscosity = 0.5; // float
        initialDensity = 1.0; // float
        force = new THREE.Vector2(0, 0); // force to apply to the simulation
    }; // class Parameters

    /**
     * Represents a single site of the lattice
     */
    static Site = class {

        // distributions Ni for 0 <= i < 9
        n = [
            0,  // i = 0
            0,  // i = 1
            0,  // i = 2
            0,  // i = 3
            0,  // i = 4
            0,  // i = 5
            0,  // i = 6
            0,  // i = 7
            0   // i = 8
        ];

        // local density (sum over Ni)
        density = 0; // float

        // local velocity u
        velocity = new THREE.Vector2(0, 0);

        
        /**
         * Default constructor, initialises the sites with initial density
         * 
         * @param {number} density Initial density at the site
         */
        constructor(density) {
            let v0 = new THREE.Vector2(0, 0);
            for(let i = 0; i < 9; ++i) {
                this.n[i] = D2Q9.NiEq(i, density, v0);
            }
        }

        /**
         * Collision step; updates density and velocity based on Ni values
         * Updates ρ and u then Ni vals
         * 
         * @param {float} omega The relaxation time inverse omega used for the simulation
         * @param {THREE.Vector2} force Forcing to add to the particles
         */
        collision(omega, force) {

            // compute ρ and u
            this.density = 0;
            this.velocity.set(0, 0);
            for(let i = 0; i < 9; ++i) {
                this.density += this.n[i];
                this.velocity.addScaledVector(D2Q9.C[i], this.n[i]);
            }
            this.velocity.multiplyScalar(1.0 / this.density); // rho * u = sum_i( c_i N_i )

            // collision step for the site
            let F = force.clone();
            F.multiplyScalar(this.density);
            for(let i = 0; i < 9; ++i) {
                let Nieq = D2Q9.NiEq(i, this.density, this.velocity);
                // Ni = Ni - omega * (Ni - Ni^eq)
                // i.e. Ni = lerp(Ni, Ni^eq, omega)
                this.n[i] = (1.0 - omega) * this.n[i] + omega * Nieq;
                // add forcing term
                // 3 Wi Ci • ρF
                this.n[i] += 3.0 * D2Q9.W[i] * D2Q9.C[i].dot(F);
            }
        }

    }; // class Site



    /// Size of the simulation domain
    width; // int
    height; // int

    /// Simulation parameters
    omega; // float
    force; // THREE.Vector2

    /// Array of individual lattice sites
    /// Front and back buffer, with ability to swap which is being read from and which is being written to
    sites; // [ D2Q9.Site[], D2Q9.Site[] ]
    flipBuffers () { this.sites = [ this.sites[1], this.sites[0] ]; } // flip front and back site buffer


    
    /**
     * Default constructor
     * 
     * @param {D2Q9.Parameters} params Simulation parameters
     */
    constructor(params) {
        this.width = params.width;
        this.height = params.height;
        this.omega = 1.0 / (3.0 * params.kinematicViscosity + 0.5); // to fulfill Navier-Stokes eq requirements (ν = ⅓(ω⁻¹ - ½))
        this.force = params.force;

        // init sites - front and back buffer
        this.sites = [[], []];
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {
                let frontSite = new D2Q9.Site(params.initialDensity);
                let backSite = new D2Q9.Site(0);
                this.sites[1].push(backSite); // front buffer for the first simulation step
                this.sites[0].push(frontSite); // back buffer - uninitialized for now
            }
        }

    }

    /**
     * Site getter, for easier access than 1d array index, accessing either a site from the front or back buffer
     * 
     * @param {number} x Index of the site to obtain on the x axis
     * @param {number} y Index of the site to obtain on the y axis
     * @param {bool} read Whether to access the site from the front or back buffer (read = front, write = back)
     */
    getSite (x, y, read) {
        while(x < 0) x += this.width; // handle negative indices
        while(y < 0) y += this.height;
        x = x % this.width; // remap to 0..width
        y = y % this.height; // remap to 0..height
        return this.sites[read ? 0 : 1][x * this.height + y];
    }

    /**
     * Single simulation timestep
     */
    step () {

        // Propagation step
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {

                // propagate from each of the 9 direction from neighbouring cells
                let site = this.getSite(x, y, false);
                for(let i = 0; i < 9; ++i) {
                    let c = D2Q9.C[i];
                    let dx = -parseInt(c.x);
                    let dy = -parseInt(c.y);
                    site.n[i] = this.getSite(x + dx, y + dy, true).n[i];
                }

            }
        }

        // Bounce-back for no-slip
        // @todo

        // Swap front and back buffers
        this.flipBuffers();

        // Collision step
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {

                this.getSite(x, y, true).collision(this.omega, this.force);

            }
        }

    }

}; // class D2Q9
