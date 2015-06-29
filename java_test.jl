type Surface

	# parameterization
	span, chord, tc, twist, sweep, dihedral, pos, nPanels

	# internal representation
	QCx,QCy,QCz # quarter chord
	LEx, TEx # leading and trailing edge
	CPx,CPy,CPz,c,t,theta,phi,ds,Lambda # control points
	eta # length along wing trace
	N # total number of panels
	symmetric

	# constants
	double rho = 1.0
	double U = 1.0
	q = 0.5*rho*U*U

	setupSurface()
  function setupSurface()

		# ----- initialize variables -----------
		N = sum(nPanels) # total number of panels
		QCx = new double[N+1] QCy = new double[N+1] QCz = new double[N+1]
		LEx = new double[N+1]
		TEx = new double[N+1]
		CPx = new double[N] CPy = new double[N] CPz = new double[N]
		c = new double[N] t = new double[N] theta = new double[N]
		phi = new double[N] ds = new double[N] Lambda = new double[N]

		int i,j,k
		int nSec = nPanels.length
		# ----------------------------------------

		# ---- initialize positions -------
		QCx[0] = pos[0] + chord[0]/4.0
		QCy[0] = pos[1]
		QCz[0] = pos[2]
		LEx[0] = pos[0]
		TEx[0] = pos[0] + chord[0]
		# --------------------------------

		# -------- compute surface variables ----------
		k = 0
		for i = 0:nSec
      dSpan = span[i]/nPanels[i]

			for j = 1:nPanels[i]

				# quarter chord locations
				QCx[k+1] = QCx[k] + dSpan*Math.tan(sweep[i])
				QCy[k+1] = QCy[k] + dSpan*Math.cos(dihedral[i])
				QCz[k+1] = QCz[k] + dSpan*Math.sin(dihedral[i])

				# leading and trailing edge locations
				double cLocal = chord[i] + j/nPanels[i]*(chord[i+1]-chord[i])
				LEx[k+1] = QCx[k+1] - cLocal/4.0
				TEx[k+1] = QCx[k+1] + 3.0/4.0*cLocal

				# control point variables
				c[k] = chord[i] + (j-0.5)/nPanels[i]*(chord[i+1]-chord[i])
				t[k] = tc[i]*chord[i] + (j-0.5)/nPanels[i]*(tc[i+1]*chord[i+1] - tc[i]*chord[i])
				theta[k] = twist[i] + (j-0.5)/nPanels[i]*(twist[i+1]-twist[i])
				phi[k] = dihedral[i]
				ds[k] = dSpan
				Lambda[k] = sweep[i]

				# control point locations
				CPx[k] = 0.5*(QCx[k] + QCx[k+1]) + 0.5*c[k]
				CPy[k] = 0.5*(QCy[k] + QCy[k+1])
				CPz[k] = 0.5*(QCz[k] + QCz[k+1])

				k++
        end
      end
		# ------------------------------------------------

		if !symmetric
			reflectXZ()
    end

    end

# ------------------------
#   Setters
# -----------------------

	function changeLength(double[] length)
		this.span = length
		setupSurface()
  end

	function changeChord(double[] chord)
		this.chord = chord
		setupSurface()
  end

	function changeTC(double[] tc)
		this.tc = tc
		setupSurface()
  end

	function changeTwist(double[] twist)
		this.twist = degToRad(twist)
		setupSurface()
  end

  function changeSweep(double[] sweep)
		this.sweep = degToRad(sweep)
		setupSurface()
  end

	function changeDihedral(double[] dihedral)
		this.dihedral = degToRad(dihedral)
		setupSurface()
  end

	function changeNumPanels(int[] nPanels)
		this.nPanels = nPanels
		setupSurface()
  end

	/**
	 * Move surface to a new position (relative to origin)
	 * @param newPosition - array with length 3 (x,y,z locations)
	 */
	function changePosition(double[] newPosition){

		double dx = (newPosition[0]-pos[0])
		double dy = (newPosition[1]-pos[1])
		double dz = (newPosition[2]-pos[2])

		add(QCx,dx) add(QCy,dy) add(QCz,dz)
		add(LEx,dx)
		add(TEx,dx)
		add(CPx,dx) add(CPy,dy) add(CPz,dz)

		pos = copy(newPosition)
	}


/* ------------------------
   Public Methods - Geometry Dependent Terms
 * ------------------------ */

	/**
	 * @return number of panels
	 */
	public int getNPanels(){
		return N
	}

	/**
	 * computes lift influence coefficient vector
	 * L = LIC'*gamma
	 * @return
	 */
	public Vector getLIC(){
		double[] LIC = new double[N]

		double factor = 1.0
		if (symmetric) factor = 2.0

		for(int i = 0 i < N i++){
			LIC[i] = factor*rho*U*Math.cos(phi[i])*ds[i]
		}

		return new Vector(LIC)
	}

	/**
	 * influence from surface s on control points on this object
	 * @param s
	 * @return induced drag influence coefficient matrix
	 */
	public Matrix getDIC(Surface s){

		int M = s.getNPanels()
		double[][] DIC = new double[N][M]
		double ny,nz,ry,rz,r2
		double factor = 1.0
		if (symmetric) factor = 2.0

		for(int i = 0 i < N i++){
			for(int j = 0 j < M j++){
				ny = QCz[i+1] - QCz[i]
				nz = QCy[i+1] - QCy[i]
				ry = s.QCy[j] - CPy[i]
				rz = s.QCz[j] - CPz[i]
				r2 = ry*ry + rz*rz
				DIC[i][j] = -factor*rho/(4.0*Math.PI*r2)*(rz*ny+ry*nz)

				# add other side of panel
				ry = s.QCy[j+1] - CPy[i]
				rz = s.QCz[j+1] - CPz[i]
				r2 = ry*ry + rz*rz
				DIC[i][j] = DIC[i][j] + factor*rho/(4.0*Math.PI*r2)*(rz*ny+ry*nz)
			}
		}

		Matrix DICmatrix = new Matrix(DIC)

		if (s.symmetric){
			Surface left = new Surface(s)
			left.flipXZroot()
			left.symmetric = false
			Matrix DICmatrix2 = this.getDIC(left)
			DICmatrix = DICmatrix.subtract(DICmatrix2) # subtract because vortices go other way on left half
		}

		return DICmatrix
	}

	/**
	 * influence from this surface on this surface
	 * Di = gamma'*DIC*gamma
	 * @return induced drag (coefficient) influence coefficient matrix
	 */
	public Matrix getDIC(){
		return getDIC(this)
	}

	/**
	 * parasite drag = cd0*q*Sref + D1'*gamma + D2'*gamma^2
	 * @param cd1
	 * @return
	 */
	public Vector getD1(double cd1){
		double factor = 1.0
		if (symmetric) factor = 2.0

		Vector D1 = new Vector(N)
		for (int i = 0 i < N i++){
			D1.set(i, factor*cd1*rho*U*ds[i])
		}
		return D1
	}

	/**
	 * parasite drag = cd0*q*Sref + D1'*gamma + D2'*gamma^2
	 * @param cd2
	 * @return
	 */
	public Vector getD2(double cd2){
		double factor = 1.0
		if (symmetric) factor = 2.0

		Vector D2 = new Vector(N)
		for (int i = 0 i < N i++){
			D2.set(i, factor*cd2*2.0*rho/c[i]/ds[i])
		}
		return D2
	}

	/**
	 * Get parasite drag using PASS method - Reynolds, Mach number dependent
	 * @param alt
	 * @param mach
	 * @param xt - transition location from 0 to 1
	 * @return
	 */
	public double getPDrag(double alt, double mach, double xt){
		double CDpS = 0

		double factor = 1.0
		if (symmetric) factor = 2.0

		for (int i = 0 i < span.length i++){
			double cr = chord[i]
			double ct = chord[i+1]
			double cbar = 0.5*(cr + ct)
			double tcbar = (tc[i]*chord[i] + tc[i+1]*chord[i+1])/(chord[i]+chord[i+1])
			double area = cbar*span[i]
			double mac = 2.0/3.0*(cr + ct - cr*ct/(cr+ct))
			double cdp = getCDpPASS(alt, mach, xt, mac, sweep[i], tcbar)
			CDpS += cdp*area
		}
		return factor*CDpS*q
	}


	/**
	 * get side force influence coefficient vector
	 * @return
	 */
	public Vector getYIC(){
		if (symmetric){
			return new Vector(N)
		} else{
			return new Vector(multiply(-rho*U,dotMultiply(sin(phi), ds)))
		}
	}


	/**
	 * gets rolling moment coefficient vector
	 * rolling moment = ROLLIC'*gamma
	 * @param cgY - y location of center of gravity
	 * @param cgZ - z location of center of gravity
	 * @return rolling moment (coefficient) influence coefficient vector
	 */
	public Vector getROLLIC(double cgY, double cgZ){
		double[] ROLLIC = new double[N]
		double factor = 1.0
		if (symmetric) factor = 0.0

		for(int i = 0 i < N i++){
			ROLLIC[i] = -factor*rho*U*ds[i]*((CPy[i]-cgY)*Math.cos(phi[i]) + (CPz[i]-cgZ)*Math.sin(phi[i]))
		}

		return new Vector(ROLLIC)
	}

	/**
	 * gets rolling moment coefficient vector about center of wing
	 * rolling moment = ROLLIC'*gamma
	 * @return rolling moment (coefficient) influence coefficient vector
	 */
	public Vector getROLLIC(){
		# assume c.g. in center of wing
		double cgY = pos[1]
		double cgZ = pos[2]
		return getROLLIC(cgY,cgZ)
	}

	/**
	 * computes pitching moment coefficient vector
	 * M = MIC'*gamma
	 * @param cgX - x location of center of gravity
	 * @return
	 */
	public Vector getMIC(double cgX){
		double[] MIC = new double[N]
		double factor = 1.0
		if (symmetric) factor = 2.0

		for(int i = 0 i < N i++){
			double qcXi = CPx[i] - 0.5*c[i]
			MIC[i] = -factor*(qcXi-cgX)*rho*U*Math.cos(phi[i])*ds[i]
		}
		return new Vector(MIC)
	}

	/**
	 * computes pitching moment coefficient vector
	 * M = MIC'*gamma
	 * assumes cg is located at the quarter chord of the surface
	 * @return
	 */
	public Vector getMIC(){
		double cgX = pos[0] + 0.25*chord[0]
		return getMIC(cgX)
	}

	/**
	 * cl = CLIC.*gamma (note: dot times)
	 * computes influence coefficients for computing local lift coefficient
	 * @return
	 */
	public Vector getCLIC(){
		double[] CLIC = new double[N]
		for(int i = 0 i < N i++){
			CLIC[i] = 2.0/U/c[i]
		}
		return new Vector(CLIC)
	}

	/**
	 * Vn = AIC*gamma
	 * get aerodynamic influence coefficient matrix from surface s
	 * on control points on this object
	 * @return
	 */
	public Matrix getAIC(Surface s){
		int M = s.getNPanels()
		double[][] u = new double[N][M]
		double[][] v = new double[N][M]
		double[][] w = new double[N][M]

		# bound vortex (BC)
		addVortexFilament(u,v,w, s.QCx,s.QCy,s.QCz,0, s.QCx,s.QCy,s.QCz,1, CPx, CPy, CPz)

		# left bound vortex (AB)
		addVortexFilament(u,v,w, s.TEx,s.QCy,s.QCz,0, s.QCx,s.QCy,s.QCz,0, CPx, CPy, CPz)

		# right bound vortex (CD)
		addVortexFilament(u,v,w, s.QCx,s.QCy,s.QCz,1, s.TEx,s.QCy,s.QCz,1, CPx, CPy, CPz)

		# trailing vortex
		double[] FFx = new double[M+1]
		for(int i = 0 i < M+1 i++){
			FFx[i] = s.TEx[i] + 1.0e9
		}
		addVortexFilament(u,v,w, FFx,s.QCy,s.QCz,0, s.TEx,s.QCy,s.QCz,0, CPx, CPy, CPz)

		# right trailing vortex
		addVortexFilament(u,v,w, s.TEx,s.QCy,s.QCz,1, FFx,s.QCy,s.QCz,1, CPx, CPy, CPz)

		# u is actually AIC - just using it because we don't actually need u (saving memory)
		for(int i = 0 i < N i++){
			for (int j = 0 j < M j++){
				u[i][j] = -v[i][j]*Math.sin(phi[i]) + w[i][j]*Math.cos(phi[i])
			}
		}

		Matrix AIC = new Matrix(u)

		if (s.symmetric){
			Surface left = new Surface(s)
			left.flipXZroot()
			left.symmetric = false
			Matrix AIC2 = this.getAIC(left)
			AIC = AIC.subtract(AIC2) #subtract because vortices opposite direction on left half
		}

		return AIC
	}

	/**
	 * Vn = AIC*gamma
	 * get aerodynamic influence coefficient matrix
	 * @return
	 */
	public Matrix getAIC(){
		return getAIC(this)
	}

	/**
	 * Gets normalwash along the surface from given freestream condition
	 * @param alpha - angle of attack in degrees
	 * @param beta - sideslip angle in degrees
	 * @return Vn - normalwash
	 */
	public Vector getVn(double alpha, double beta){
		alpha = alpha/180.0*Math.PI
		beta = beta/180.0*Math.PI

		double[] Vn = new double[N]
		for(int i = 0 i < N i++){
			Vn[i] = U*(Math.cos(alpha)*Math.cos(beta)*Math.sin(theta[i])
					   + Math.sin(beta)*Math.cos(theta[i])*Math.sin(phi[i])
					   + Math.sin(alpha)*Math.cos(beta)*Math.cos(theta[i])*Math.cos(phi[i]))
		}
		return new Vector(Vn)
	}

	/**
	 * get just the geometry dependent terms of Vn
	 * Vn = Vn1*cos(alpha)*cos(beta) + Vn2*sin(beta) + Vn3*sin(alpha)*cos(beta)
	 * @param Vn1
	 * @param Vn2
	 * @param Vn3
	 */
	public void getVnComponents(Vector Vn1, Vector Vn2, Vector Vn3){

		for(int i = 0 i < N i++){
			Vn1.set(i, U*Math.sin(theta[i]))
			Vn2.set(i, U*Math.cos(theta[i])*Math.sin(phi[i]))
			Vn3.set(i, U*Math.cos(theta[i])*Math.cos(phi[i]))
		}
		return
	}

	/**
	 * Get derivative of Vn w.r.t theta (twist)
	 */
	public Vector getdVndTheta(double alpha, double beta){
		alpha = alpha/180.0*Math.PI
		beta = beta/180.0*Math.PI

		double[] dVn = new double[N]
		for(int i = 0 i < N i++){
			dVn[i] = U*(Math.cos(alpha)*Math.cos(beta)*Math.cos(theta[i])
					   - Math.sin(beta)*Math.sin(theta[i])*Math.sin(phi[i])
					   - Math.sin(alpha)*Math.cos(beta)*Math.sin(theta[i])*Math.cos(phi[i]))
		}
		return new Vector(dVn)
	}


	/**
	 * compute external normalwash from an external velocity field
	 * @param field
	 * @return
	 */
	public Vector getVnExternal(VelocityField field){
		double[] VnExternal = new double[N]
		for(int i = 0 i < N i++){
			VnExternal[i] = field.getVx(CPx[i], CPy[i], CPz[i])*Math.sin(theta[i])
						  - field.getVy(CPx[i], CPy[i], CPz[i])*Math.cos(theta[i])*Math.sin(phi[i])
						  + field.getVz(CPx[i], CPy[i], CPz[i])*Math.cos(theta[i])*Math.cos(phi[i])
		}
		return new Vector(VnExternal)
	}


	/**
	 * compute induced drag due to normalwash from an external velocity field
	 * @param factor - multiplies number of control points by factor for more accurate integration
	 * @param field
	 * @return
	 */
	public double getDiExternal(int factor, VelocityField field, double[] gamma, int position){

		int nInterp = factor*N

		computeEta()
		double eta_end = eta[N-1] + ds[N-1]/2.0
		double[] etaInterp = linspace(0.0,eta_end,nInterp)
		double ds = eta_end/(nInterp-1)

		# TODO: there is a better way to do this.  why interpolate.  I know what it should be exactly.
		# interpolate
		double[] x = interp1(eta,CPx,etaInterp)
		double[] y = interp1(eta,CPy,etaInterp)
		double[] z = interp1(eta,CPz,etaInterp)
		double[] twist = interp1(eta,theta,etaInterp)
		double[] dihedral = interp1(eta,phi,etaInterp)
		double[] g = interp1(eta,gamma,etaInterp)

		if (position == 0) g[0] = 0
		else if (position == 1) g[nInterp-1] = 0

		double VnInterp
		double Di = 0

		# compute new influence coefficient vector
		for (int i = 0 i < nInterp i++){
			VnInterp = field.getVx(x[i],y[i],z[i])*Math.sin(twist[i])
     			     - field.getVy(x[i],y[i],z[i])*Math.cos(twist[i])*Math.sin(dihedral[i])
     				 + field.getVz(x[i],y[i],z[i])*Math.cos(twist[i])*Math.cos(dihedral[i])
			Di += -rho*VnInterp*g[i]*ds
		}

		return Di
	}




	/**
	 * @return get interpolated set of y values at factor times the current panel density
	 */
	public Vector getYinterp(int factor){
		computeEta()
		double eta_end = eta[N-1] + ds[N-1]/2.0
		double[] etaInterp = linspace(0.0,eta_end,N*factor)
		return new Vector(interp1(eta,CPy,etaInterp))
	}

	/**
	 * @return get interpolated set of y values at factor times the current panel density
	 */
	public Vector getZinterp(int factor){
		computeEta()
		double eta_end = eta[N-1] + ds[N-1]/2.0
		double[] etaInterp = linspace(0.0,eta_end,N*factor)
		return new Vector(interp1(eta,CPz,etaInterp))
	}

	public Vector getGammaInterp(int factor, double[] gamma){
		computeEta()
		double eta_end = eta[N-1] + ds[N-1]/2.0
		double[] etaInterp = linspace(0.0,eta_end,N*factor)
		return new Vector(interp1(eta,gamma,etaInterp))
	}


	/**
	 * deflect ailerons antisymmetrically
	 * @param aileron - angle to deflect.
	 * @param start - start of the control surface in % semi-span starting from root
	 * @param finish - end of the control surface in % semi-span
	 * @return
	 */
	public void deflectAilerons(double aileron, double start, double finish){
		if (symmetric){
			System.out.println("Symmetric aircraft is already trimmed")
			return
		}

		aileron = aileron*Math.PI/180.0

		int n1 = (int) (N/2.0*(1.0-finish))
		int n2 = (int) (n1 + N/2.0*(finish-start))
		int n3 = (int) (n2 + N*start)
		int n4 = (int) (n3 + N/2.0*(finish-start))
		int i

		for (i = n1 i < n2 i++){
			theta[i] = theta[i] + aileron
		}
		for (i = n3 i < n4 i++){
			theta[i] = theta[i] - aileron
		}
	}



	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Surface [chord=" + Arrays.toString(chord) + ", dihedral="
				+ Arrays.toString(dihedral) + ", nPanels="
				+ Arrays.toString(nPanels) + ", pos=" + Arrays.toString(pos)
				+ ", span=" + Arrays.toString(span) + ", sweep="
				+ Arrays.toString(sweep) + ", symmetric=" + symmetric + ", tc="
				+ Arrays.toString(tc) + ", theta=" + Arrays.toString(theta)
				+ ", twist=" + Arrays.toString(twist) + "]"
	}

	/* ------------------------
   Private Helper Methods
 * ------------------------ */
	/*
	 * convert array of values in degrees to radians
	 */
	private static double[] degToRad(double[] deg){
		double[] rad = new double[deg.length]

		for (int i = 0 i < deg.length i++){
			rad[i] = deg[i]/180.0*Math.PI
		}
		return rad
	}

	/*
	 * sum up values in integer array
	 */
	private static int sum(int[] a){
		int sum = 0
		for (int i: a) sum += i

		return sum
	}

	/*
	 * Creates a full wing rather than just a half wing.  (reflects about the
	 * X-Z plane).
	 */
	private void reflectXZ(){

		Surface left = new Surface(this)
		left.flipXZtip()

		QCx = combineMerge(left.QCx,QCx)
		QCy = combineMerge(left.QCy,QCy)
		QCz = combineMerge(left.QCz,QCz)

		LEx = combineMerge(left.LEx,LEx)
		TEx = combineMerge(left.TEx,TEx)

		CPx = combine(left.CPx,CPx)
		CPy = combine(left.CPy,CPy)
		CPz = combine(left.CPz,CPz)
		c = combine(left.c,c)
		t = combine(left.t,t)
		theta = combine(left.theta,theta)
		phi = combine(left.phi,phi)
		ds = combine(left.ds,ds)
		Lambda = combine(left.Lambda,Lambda)

		N = 2*N
	}

	/*
	 * Flip wing about the X-Z plane. have 0 start from root
	 */
	private void flipXZroot(){

		# first move wing to y = 0
		double[] savePosition = copy(pos)
		double[] newPosition = {pos[0],0,pos[2]}
		changePosition(newPosition)

		QCy = negative(QCy)
		CPy = negative(CPy)
		phi = negative(phi)

		# move back to where it started
		changePosition(savePosition)
	}


	/*
	 * Flip wing about the X-Z plane. have 0 start from tip
	 */
	private void flipXZtip(){

		# first move wing to y = 0
		double[] savePosition = copy(pos)
		double[] newPosition = {pos[0],0,pos[2]}
		changePosition(newPosition)

		# -------- quarter chord, leading and trailing edge ---------
		QCx = flip(QCx)
		QCy = negative(flip(QCy))
		QCz = flip(QCz)
		LEx = flip(LEx)
		TEx = flip(TEx)
		# --------------------------------------------------

		# ---------- control points -----------------
		CPx = flip(CPx)
		CPy = negative(flip(CPy))
		CPz = flip(CPz)
		c = flip(c)
		t = flip(t)
		theta = flip(theta)
		phi = negative(flip(phi))
		ds = flip(ds)
		Lambda = flip(Lambda)
		# ----------------------------------------

		# move back to where it started
		changePosition(savePosition)
	}

	/*
	 * combine two arrays into one large one
	 * [d1 d2]
	 * @param d1
	 * @param d2
	 * @return
	 */
	private double [] combine(double[] d1, double[] d2){
		int n1 = d1.length
		int n2 = d2.length
		double[] dnew = new double[n1 + n2]
		System.arraycopy(d1, 0, dnew, 0, n1)
		System.arraycopy(d2, 0, dnew, n1, n2)
		return dnew
	}

	/*
	 * combine two arrays but merge the middle part (assumed they are equal)
	 * [d1 d2(2:end)]
	 * @param d1
	 * @param d2
	 * @return
	 */
	private double[] combineMerge(double[] d1, double[] d2){
		int n1 = d1.length
		int n2 = d2.length
		double[] dnew = new double[n1 + n2 -1]
		System.arraycopy(d1, 0, dnew, 0, n1)
		System.arraycopy(d2, 1, dnew, n1, n2-1)
		return dnew
	}

	/*
	 * computes velocity from a vortex filaments starting at rA, ending at rB at control
	 * points rC.  velocity is added to whatever is passed into v and w.
	 *
	 * u,v,w(i,j) = axial,lateral and vertical velocity induced at control point rC(i)
	 * due to vortex filament that goes from rA(j) to rB(j)
	 *
	 * see Bertin & Smith
	 *
	 */
	private void addVortexFilament(double[][]u, double[][]v, double[][]w,
			double[] xA, double[] yA, double[] zA, int startA,
			double[] xB, double[] yB, double[] zB, int startB,
			double[] xC, double[] yC, double[] zC){

		int m = yC.length
		int n = yA.length-1
		double x1,x2,x21,y1,y2,y21,z1,z2,z21,denom1,i1,j1,k1,frac1_1,frac1_2,frac1,factor
		int p, q

		for (int i = 0 i < m  i++){
			p = startA q = startB
			for (int j = 0 j < n j++){
				x1 = xC[i] - xA[p]
		        x2 = xC[i] - xB[q]
		        x21 = xB[q] - xA[p]
		        y1 = yC[i] - yA[p]
		        y2 = yC[i] - yB[q]
		        y21 = yB[q] - yA[p]
		        z1 = zC[i] - zA[p]
		        z2 = zC[i] - zB[q]
		        z21 = zB[q] - zA[p]
		        p++ q++

		        denom1 = Math.pow(y1*z2-y2*z1,2) + Math.pow(x1*z2-x2*z1,2) + Math.pow(x1*y2-x2*y1,2)
		        i1 = y1*z2 - y2*z1
		        j1 = -x1*z2 + x2*z1
		        k1 = x1*y2 - x2*y1
		        frac1_1 = (x21*x1 + y21*y1 + z21*z1)/Math.sqrt(x1*x1 + y1*y1 + z1*z1)
		        frac1_2 = (x21*x2 + y21*y2 + z21*z2)/Math.sqrt(x2*x2 + y2*y2 + z2*z2)
		        frac1 = frac1_1 - frac1_2
		        factor = 1.0/(4.0*Math.PI)/denom1*frac1

		        u[i][j] = u[i][j] + i1*factor
		        v[i][j] = v[i][j] + j1*factor
		        w[i][j] = w[i][j] + k1*factor
			}
		}

		return
	}

	private void computeEta(){

		eta = new double[N]
		eta[0] = ds[0]/2.0
		for(int i = 1 i < N i++){
			eta[i] = eta[i-1] + ds[i-1]/2.0 + ds[i]/2.0
		}
	}

	private double getCDpPASS(double alt, double mach, double xt, double mac, double sweep, double tc){

		# skin friction drag coefficient
		double Cf = Atmosphere.cf(alt, mach, mac, xt)

		# form factor
		double cossw = Math.cos(sweep)
		double m0 = 0.5
		double z = (2-m0*m0)*cossw/Math.sqrt(1-Math.pow(m0*cossw,2))
		double k = 1 + tc*z + Math.pow(tc,4)*100.0

		# wetted area / S
		double SwetS = 2.0*(1+0.2*tc)

		# parasite drag
		return Cf*k*SwetS
	}

#	private double getCDcPASS(double CL, double sweep, double tc, double mach, double supercrit){
#		  double cossw = Math.cos(sweep)
#	      double clp   = Math.abs(CL/ (cossw*cossw))
#	      double tcp   = tc / cossw
#
#	      double mcc = .954-.235*clp+.0259*clp*clp
#	      mcc = mcc-(1.963-1.078*clp+.350*clp*clp)*tcp
#	      mcc = mcc+(2.969-2.738*clp+1.469*clp*clp)*tcp*tcp
#
#	      # correction for scw
#	      mcc = mcc + supercrit*.06
#	      mcc = mcc / cossw
#
#	      double rm = mach / mcc
#	      double dm = rm-1.
#
#	      double cdc = 0.
#	      if (rm >= .5 && rm < .8) cdc = 1.3888952e-4+5.5556e-4*dm+5.5556192e-4*dm*dm
#	      if (rm >= .8 && rm < .95) cdc = 7.09e-4+6.733e-3*dm+.01956*dm*dm+.01185*dm*dm*dm
#	      if (rm >= .95 && rm < 1.0) cdc = .001000+.02727*dm+.4920*dm*dm+3.57385*dm*dm*dm
#	      if (rm >= 1.0) cdc = .001000+.02727*dm-.1952*dm*dm+19.09*dm*dm*dm
#
#	      return cdc * Math.pow(cossw,3)
#	}

#	/*
#	 * convert string to double array
#	 */
#	private static double[] stringToDoubles(String string){
#		String[] stringArray = string.split("\\s+")
#		int n = stringArray.length
#
#		double[] result = new double[n]
#		for (int i = 0 i < n i++){
#			result[i] = Double.parseDouble(stringArray[i])
#		}
#
#		return result
#	}
#
#	/*
#	 * convert string to integer array
#	 */
#	private static int[] stringToInts(String string){
#		String[] stringArray = string.split("\\s+")
#		int n = stringArray.length
#
#		int[] result = new int[n]
#		for (int i = 0 i < n i++){
#			result[i] = Integer.parseInt(stringArray[i])
#		}
#
#		return result
#	}









	# ---------------------------------------------------------



}
