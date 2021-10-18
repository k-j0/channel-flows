
/**
 * D2Q9 lattice boltzmann method simulation on rectangular lattice
 */
class D2Q9 {

    /// Simulation constants (always the same across any D2Q9 sim)
    // Order of i is chosen clockwise, to make full bounce-back conditions easy to compute with modulos
    // 6 7 8
    // 5 0 1
    // 4 3 2
    // Hence the weights and velocities look 'out of order' compared to typical order with 1 <= i <= 4 on the sides and 5 <= i <= 8 on the corners
    static W = [ // weights
        4.0 / 9.0,  // i = 0 (centre)
        1.0 / 9.0,  // i = 1 (side)
        1.0 / 36.0, // i = 2 (corner)
        1.0 / 9.0,  // i = 3 (side)
        1.0 / 36.0, // i = 4 (corner)
        1.0 / 9.0,  // i = 5 (side)
        1.0 / 36.0, // i = 6 (corner)
        1.0 / 9.0,  // i = 7 (side)
        1.0 / 36.0  // i = 8 (corner)
    ];
    static C = [ // velocities
        new THREE.Vector2(0, 0),    // i = 0
        new THREE.Vector2(1, 0),    // i = 1
        new THREE.Vector2(1, 1),    // i = 2
        new THREE.Vector2(0, 1),    // i = 3
        new THREE.Vector2(-1, 1),   // i = 4
        new THREE.Vector2(-1, 0),   // i = 5
        new THREE.Vector2(-1, -1),  // i = 6
        new THREE.Vector2(0, -1),   // i = 7
        new THREE.Vector2(1, -1)    // i = 8
    ];
    

    /**
     * Computes the equilibrium value Ni^eq for a site
     * 
     * @param {number} i Direction for which to get Ni^eq
     * @param {number} density The local density ρ
     * @param {THREE.Vector2} velocity The velocity u
     */
    static NiEq (i, density, velocity) {
        let uDotU = velocity.dot(velocity); // u • u
        let ciDotU = D2Q9.C[i].dot(velocity); // c_i • u
        return D2Q9.W[i] * density * (1.0 + 3.0 * ciDotU + 4.5 * ciDotU * ciDotU - 1.5 * uDotU);
    }

    /**
     * Computes the opposite direction from direction i
     * 
     * @param {number} i Direction for which to get the opposing dir
     * @returns {number} Opposite direction to input
     */
    static opposite (i) {
        if (i == 0) return 0;
        return ((i - 1 + 4) % 8) + 1;
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

        // whether the site is solid (boundary site); particles that end up on solid sites will bounce back
        boundary = false;

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

        // temp distributions Ni after collision step
        nt = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0
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
         * Updates ρ and u then Ni vals (nt)
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
                this.velocity.addScaledVector(D2Q9.C[i], this.n[i]); // ρu = sum_i( c_i * N_i )
            }
            this.velocity.x += 1e-2; // horizontal pressure gradient; TODO remove this
            this.velocity.multiplyScalar(1.0 / this.density); // ρu = sum_i( c_i * N_i ) hence u = sum_i(c_i * N_i ) / ρ

            // collision step for the site
            let F = force.clone();
            F.multiplyScalar(this.density);
            for(let i = 0; i < 9; ++i) {
                let Nieq = D2Q9.NiEq(i, this.density, this.velocity);
                // Ni = Ni - omega * (Ni - Ni^eq)
                // i.e. Ni = lerp(Ni, Ni^eq, omega)
                this.nt[i] = (1.0 - omega) * this.n[i] + omega * Nieq;
                if (this.nt[i] < 0) {
                    this.nt[i] = 0;
                }
                // add forcing term   3 Wi Ci • F   with   F = ρf
                //this.n[i] += 3.0 * D2Q9.W[i] * D2Q9.C[i].dot(F); // <- TODO add this back in!
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
    sites; // D2Q9.Site[]
    
    /// Time-dependent values (not needed for simulation, just for display purposes)
    t = 0; // time, i.e. number of steps taken
    density = 0; // average density over the lattice
    velocity = new THREE.Vector2(0, 0); // average velocity over the lattice


    
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

        // init sites
        this.sites = [];
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {
                this.sites.push(new D2Q9.Site(params.initialDensity));
            }
        }

    }

    /**
     * Adds a boundary on a single site
     * 
     * @param {int} x X-coordinate of the site to make into a boundary site
     * @param {int} y Y-coordinate of the site to make into a boundary site
     */
    addBoundary (x, y) {
        this.getSite(x, y).boundary = true;
        for (let i = 0; i < 9; ++i) {
            this.getSite(x, y).n[i] = 0;
        }
    }

    /**
     * Adds a horizontal bounce-back boundary for sites at y
     * 
     * @param {int} y Y-coordinate of sites to make into boundaries
     */
    addHorizontalBoundary(y) {
        for(let x = 0; x < this.width; ++x) {
            this.addBoundary(x, y);
        }
    }

    /**
     * Site getter, for easier access than 1d array index
     * Handles periodicity on both axes
     * 
     * @param {number} x Index of the site to obtain on the x axis
     * @param {number} y Index of the site to obtain on the y axis
     */
    getSite (x, y) {
        while(x < 0) x += this.width; // handle negative indices
        while(y < 0) y += this.height;
        x = x % this.width; // remap to 0..width
        y = y % this.height; // remap to 0..height
        return this.sites[x * this.height + y];
    }

    /**
     * Single simulation timestep
     */
    step () {

        this.density = 0;
        this.velocity.set(0, 0);

        // Collision step
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // get site to update
                let site = this.getSite(x, y);
                
                // skip calculations for obstacle sites
                if (site.boundary) continue;

                // collision step for the site
                site.collision(this.omega, this.force);
            }
        }

        // Propagation step
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // get site to update
                let site = this.getSite(x, y);

                // skip boundary sites
                if (site.boundary) continue;

                // stream particle distribs from neighbours
                for (let i = 0; i < 9; ++i) {
                    let dx = D2Q9.C[i].x;
                    let dy = D2Q9.C[i].y;
                    site.n[i] = this.getSite(x + dx, y + dy).nt[i];
                }

            }
        }

        // Bounce-back boundary conditions
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // get site to push particles away from
                let site = this.getSite(x, y);
                
                // skip fluid sites
                if (!site.boundary) continue;

                // half-way bounce-back
                for (let i = 0; i < 9; ++i) {
                    let dx = -D2Q9.C[i].x;
                    let dy = -D2Q9.C[i].y;
                    site.n[i] = this.getSite(x + dx, y + dy).n[D2Q9.opposite(i)];
                }
            }
        }

        // Compute avg density and velocity
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                let site = this.getSite(x, y);
                this.density += site.density;
                this.velocity.add(site.velocity);

            }
        }
        this.density /= this.width * this.height;
        this.velocity.multiplyScalar(1.0 / (this.width * this.height));

        // t -> t+1
        ++this.t;

    }

}; // class D2Q9
