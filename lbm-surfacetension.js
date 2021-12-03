
/**
 * D2Q9 lattice boltzmann method simulation on rectangular lattice with surface tension
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
     * Computes the lattice tensor Q' applied to a vector r
     * Q'_{iαβ} r_α r_β
     * 
     * @param {int} i Site direction to apply
     * @param {THREE.Vector2} r Vector to apply the lattice tensor to
     * @returns {float} Result of the tensor operation
     */
    static Q (i, r) {
        let ciDotR = D2Q9.C[i].dot(r); // c_i • r
        return 1.5 * ciDotR * ciDotR - 0.5 * r.dot(r); // 3/2 (c_i • r)^2 - 1/2 r • r
    }

    /**
     * Computes the equilibrium value Ni^eq for a site
     * 
     * @param {number} i Direction for which to get Ni^eq
     * @param {number} density The local density ρ
     * @param {THREE.Vector2} velocity The velocity u
     */
    static NiEq (i, density, velocity) {
        let ciDotU = D2Q9.C[i].dot(velocity); // c_i • u
        let nieq = D2Q9.W[i] * density * (1.0 + 3.0 * ciDotU + 3.0 * D2Q9.Q(i, velocity));
        if (nieq < 0) return 0; // would only happen due to numerical imprecisions
        return nieq;
    }

    /**
     * Computes the omega value corresponding to a kinematic viscosity to fulfill Navier-Stokes eq requirements (ν = ⅓(ω⁻¹ - ½))
     * 
     * @param {number} nu Kinematic viscosity to return omega for 
     * @returns {number} Omega value corresponding to the kinematic viscosity
     */
    static omega (nu) {
        return 1.0 / (3.0 * nu + 0.5);
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
        kinematicViscosity1 = 0.5; // float, kinematic viscosity of the first fluid
        kinematicViscosity2 = 0.5; // float, kinematic viscosity of the second fluid
        force = (time) => new THREE.Vector2(0, 0); // force to apply to the simulation; function: float time -> THREE.Vector2
        surfaceTension = 0.0; // float, free parameter A
        interfaceWidth = 1.0; // float, beta parameter to define interface width for surface tension calculation
        initialDistribution = (x, y, r, b) => {
            const density = 3.0;
            const v0 = new THREE.Vector2(0, 0);
            const d = Math.random(); // random distribution between red and blue
            for(let i = 0; i < 9; ++i) {
                // each fluid gets a part of the original Ni depending on the distribution d
                r[i] = D2Q9.NiEq(i, density, v0) * d;
                b[i] = D2Q9.NiEq(i, density, v0) * (1.0 - d);
            }
        }; // function: int x, int y, float[9] r, float[9] b
    }; // class Parameters

    /**
     * Represents a single site of the lattice
     */
    static Site = class {

        // whether the site is solid (boundary site); particles that end up on solid sites will bounce back
        boundary = false;

        // distributions Ri and Bi for 0 <= i < 9
        r = [
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
        b = [
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

        // temp distributions Ri and Bi after collision step
        rt = [
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
        bt = [
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

        // local colour gradient f
        colourGradient = new THREE.Vector2(0, 0);

        // local kinematic viscosity ν
        nu = 0; // float

        
        /**
         * Default constructor, initialises the sites with initial density and 0 velocity
         * 
         * @param {int} x The site's x coordinate
         * @param {int} y The site's y coordinate
         * @param {function} distributionFn Distribution function to use to fill the site initially
         */
        constructor(x, y, distributionFn) {
            distributionFn(x, y, this.r, this.b);
        }

        /**
         * Returns N_i for the site as sum of B_i and R_i
         * 
         * @param {int} i The direction for which to get N_i
         * @returns {float} Quantity N_i for the site
         */
        N(i) {
            return this.r[i] + this.b[i];
        }

        /**
         * Collision step; updates density and velocity based on Ni values
         * Updates ρ & u for the site, then its Ri and Bi vals (into rt/bt instead of r/b)
         * 
         * @param {float} nu1 Kinematic viscosity of the first coloured fluid
         * @param {float} nu2 Kinematic viscosity of the second coloured fluid
         * @param {THREE.Vector2} force Force to add to the particles
         * @param {float} a Free parameter used to set the surface tension
         * @param {float} beta Surface tension interface width
         */
        collision(nu1, nu2, force, a, beta) {

            // compute ρ and u
            this.density = 0;
            this.velocity.set(0, 0);
            let n = [0, 0, 0, 0, 0, 0, 0, 0, 0];
            for(let i = 0; i < 9; ++i) {
                n[i] = this.N(i);
                this.density += n[i];
                this.velocity.addScaledVector(D2Q9.C[i], n[i]); // ρu = sum_i( c_i * N_i )
            }
            this.velocity.multiplyScalar(1.0 / this.density); // ρu = sum_i( c_i * N_i ) hence u = sum_i(c_i * N_i ) / ρ

            // collision step for the site
            let F = new THREE.Vector2(force.x, force.y);
            F.x *= this.density;
            F.y *= this.density;
            let nPrime = [0, 0, 0, 0, 0, 0, 0, 0, 0];
            let R = 0, B = 0;
            for (let i = 0; i < 9; ++i) {
                R += this.r[i];
                B += this.b[i];

                // compute Ni^{eq}
                let Nieq = D2Q9.NiEq(i, this.density, this.velocity);

                // compute relaxation parameter
                if (this.r[i] == 0) this.nu = nu2; // only contribution from the blue fluid
                else if (this.b[i] == 0) this.nu = nu1; // only contribution from the red fluid
                else {
                    let t = this.r[i] / n[i]; // normalized value representing the fraction of R_i compared to the total quantity (R_i + B_i)
                    this.nu = nu1 * t + nu2 * (1.0 - t); // lerp between nu2 at t == 0 and nu1 at t == 1
                }
                let omega = D2Q9.omega(this.nu);

                // compute N'i
                // with forcing term   3 Wi Ci • F   with   F = ρf
                // and colour gradient-based forcing A |f| W_i Q'f/|f|
                let forcingTerm = 3.0 * D2Q9.W[i] * D2Q9.C[i].dot(F); // external forcing term (outside colour gradient)
                let fHat = this.colourGradient.clone().normalize(); // forcing from colour gradient
                nPrime[i] = n[i] - omega * (n[i] - Nieq) + forcingTerm + a * this.colourGradient.length() * D2Q9.W[i] * D2Q9.Q(i, fHat);

                // ensure numerical correctness even with imprecisions
                if (nPrime[i] < 0) nPrime[i] = 0;

            }

            // Antidiffusive recolouring step
            for (let i = 0; i < 9; ++i) {

                // compute Nieq(ρ, 0)
                let Nieq0 = D2Q9.NiEq(i, this.density, new THREE.Vector2(0, 0));
                let cosPhi = 0; // force cosφ = 0 for c_i = 0
                let fLen = this.colourGradient.length();
                if (i > 0 && fLen > 0) {
                    cosPhi = D2Q9.C[i].dot(this.colourGradient) / fLen;
                }
                
                this.rt[i] = R / (R + B) * nPrime[i] + beta * R * B * Nieq0 * cosPhi / ((R+B)*(R+B))
                this.bt[i] = B / (R + B) * nPrime[i] - beta * R * B * Nieq0 * cosPhi / ((R+B)*(R+B));

                // fix numerical imprecisions
                if(this.rt[i] < 0) this.rt[i] = 0;
                if(this.bt[i] < 0) this.bt[i] = 0;
            }
        }

        /**
         * Make the site into a boundary/obstacle
         */
        setBoundary() {
            this.boundary = true;
            for (let i = 0; i < 9; ++i) {
                this.r[i] = this.b[i] = 0;
            }
        }

    }; // class Site



    /// Size of the simulation domain
    width; // int
    height; // int

    /// Simulation parameters
    nu1; // float
    nu2; // float
    force; // function: float time -> THREE.Vector2
    a; // float; free parameter controlling surface tension
    beta; // float; controls surface tension interface width

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
        this.nu1 = params.kinematicViscosity1;
        this.nu2 = params.kinematicViscosity2;
        this.force = params.force;
        this.a = params.surfaceTension;
        this.beta = params.interfaceWidth;

        // init sites
        this.sites = [];
        for(let x = 0; x < this.width; ++x) {
            for(let y = 0; y < this.height; ++y) {
                this.sites.push(new D2Q9.Site(x, y, params.initialDistribution));
            }
        }

    }

    /**
     * Adds a horizontal bounce-back boundary for sites at y
     * 
     * @param {int} y Y-coordinate of sites to make into boundaries
     */
    addHorizontalBoundary(y) {
        for(let x = 0; x < this.width; ++x) {
            this.getSite(x, y).setBoundary();
        }
    }

    /**
     * Adds a vertical bounce-back boundary for sites at x
     * 
     * @param {int} x X-coordinate of sites to make into boundaries
     */
    addVerticalBoundary(x) {
        for(let y = 0; y < this.height; ++y) {
            this.getSite(x, y).setBoundary();
        }
    }

    /**
     * Adds a circular bounce-back boundary for sites centred at x, y with radius r
     * 
     * @param {int} x X-coordinate of the centre of the circle
     * @param {int} y Y-coordinate of the centre of the circle
     * @param {float} r Radius of the circle 
     */
    addCircularBoundary(x, y, r) {
        for(let u = 0; u < this.width; ++u) {
            for(let v = 0; v < this.height; ++v) {
                let dx = x - u;
                let dy = y - v;
                if (dx*dx + dy*dy <= r*r) {
                    this.getSite(u, v).setBoundary();
                }
            }
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

        // Update colour gradient at each site
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // update colour gradient at the site based on neighbours
                // f = sum_i c_i sum_j R_j(x + c_i) - B_j(x + c_i)
                let site = this.getSite(x, y);
                site.colourGradient.set(0, 0);
                for (let i = 0; i < 9; ++i) {
                    let dx = D2Q9.C[i].x;
                    let dy = D2Q9.C[i].y;
                    for (let j = 0; j < 9; ++j) {
                        let neighbour = this.getSite(x - dx, y - dy);
                        let rj = neighbour.r[j];
                        let bj = neighbour.b[j];
                        site.colourGradient.add(D2Q9.C[i].clone().multiplyScalar(rj - bj));
                    }
                }

            }
        }

        // Collision step + boundary conditions
        let f = this.force(this.t);
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // get site to update
                let site = this.getSite(x, y);

                // collision step for the site
                site.collision(this.nu1, this.nu2, f, this.a, this.beta);

                // Compute avg density and velocity across lattice
                // This isn't required, but we want to display it on screen
                if(!site.boundary) {
                    this.density += site.density;
                    this.velocity.add(site.velocity);
                }
                
                // half-way bounce-back
                // skip fluid sites for bounce-back
                if (site.boundary) {
                    for (let i = 0; i < 9; ++i) {
                        site.rt[i] = site.r[D2Q9.opposite(i)];
                        site.bt[i] = site.b[D2Q9.opposite(i)];
                    }
                }
            }
        }

        // Average out these values
        this.density /= this.width * this.height;
        this.velocity.multiplyScalar(1.0 / (this.width * this.height));

        // Propagation/streaming step
        for (let x = 0; x < this.width; ++x) {
            for (let y = 0; y < this.height; ++y) {

                // get site to update
                let site = this.getSite(x, y);

                // stream particle distribs from neighbours
                for (let i = 0; i < 9; ++i) {
                    let dx = D2Q9.C[i].x;
                    let dy = D2Q9.C[i].y;
                    site.r[i] = this.getSite(x + dx, y + dy).rt[i];
                    site.b[i] = this.getSite(x + dx, y + dy).bt[i];
                }

            }
        }

        // t -> t+1
        ++this.t;
    }





    // ---------------------------------------------------------------------------------------------




    /**
     * Returns the maximum x velocity in the channel at the current point in time
     */
    getMaxXVelocity () {
        let m = Number.NEGATIVE_INFINITY;
        for(let y = 0; y < this.height; ++y) {
            for (let x = 0; x < this.width; ++x) {
                let site = this.getSite(x, y);
                if(site.boundary) continue;
                if(site.velocity.x > m) m = site.velocity.x;
            }
        }
        return m;
    }

    /**
     * Returns the x velocity as a function of y as discrete data points across the channel
     */
    getXVelocityData () {
        let data = [];
        for (let y = 0; y < this.height; ++y) {
            // average through channel
            let velX = 0;
            for (let x = 0; x < this.width; ++x) {
                let site = this.getSite(x, y);
                if(!site.boundary) {
                    velX += site.velocity.x;
                }
            }
            velX /= this.width;
            data.push([y - (this.height - 1) * 0.5, velX]);
        }

        return data;
    }

    /**
     * Returns the surface tension free-energy between two lattice sites as σ = R/6 ∆ρ
     * 
     * @param {float} curvature Radius of curvature R
     * @param {int} x1 X-position of the first site
     * @param {int} y1 Y-position of the first site
     * @param {int} x2 X-position of the second site
     * @param {int} y2 Y-position of the second site
     */
    getSurfaceTension (curvature, x1, y1, x2, y2) {
        let deltaRho = Math.abs(this.getSite(x2, y2).density - this.getSite(x1, y1).density);
        return curvature * deltaRho / 6.0;
    }

    /**
     * Plots the x velocity as a function of y (across channel)
     * 
     * @param {string} colour Colour to use on the chart
     */
    plotXVelocity (colour = 'black') {

        let velData = this.getXVelocityData();
        let labels = velData.map((a) => a[0]);
        let data = velData.map((a) => a[1]);

        plot(labels, data, 'Velocity across channel at t = ' + this.t, colour);
    }

    /**
     * Plots the average values of N_i across the entire lattice
     * 
     * @param {string} colour Colour to use on the chart
     */
    plotOccupations (colour = 'red') {
        
        let labels = [];
        let data = [];

        for (let i = 0; i < 9; ++i) {
            // average through lattice
            let avg = 0;
            for (let x = 0; x < this.width; ++x) {
                for (let y = 0; y < this.height; ++y) {
                    let site = this.getSite(x, y);
                    if (!site.boundary) {
                        avg += site.n[i];
                    }
                }
            }
            avg /= this.width * this.height * D2Q9.W[i];
            labels.push('i = ' + i);
            data.push(avg);
        }

        plot(labels, data, 'Occupation average across lattice at t = ' + this.t, colour);
    }

    /**
     * Returns the amount of reds across the center of the lattice
     */
    getRedAmountData () {
        let data = [];
        const y = this.height / 2;
        for (let x = 0; x < this.width; ++x) {
            let r = 0;
            for (let i = 0; i < 9; ++i) {
                r += this.getSite(x, y).r[i];
            }
            data.push([x, r]);
        }
        return data;
    }

    /**
     * Plots the R(x) across the center of the lattice
     * 
     * @param {string} colour Colour to use on the chart
     */
    plotReds (colour = 'red') {

        let redData = this.getRedAmountData();
        let labels = redData.map((a) => a[0]);
        let data = redData.map((a) => a[1]);

        plot(labels, data, 'Red quantity across lattice centre at t = ' + this.t, colour);
        
    }

}; // class D2Q9
