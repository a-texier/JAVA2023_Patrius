package progmission;

import fr.cnes.sirius.patrius.assembly.models.SensorModel;
import fr.cnes.sirius.patrius.attitudes.*;
import fr.cnes.sirius.patrius.events.CodedEvent;
import fr.cnes.sirius.patrius.events.CodedEventsLogger;
import fr.cnes.sirius.patrius.events.GenericCodingEventDetector;
import fr.cnes.sirius.patrius.events.Phenomenon;
import fr.cnes.sirius.patrius.events.postprocessing.AndCriterion;
import fr.cnes.sirius.patrius.events.postprocessing.ElementTypeFilter;
import fr.cnes.sirius.patrius.events.postprocessing.Timeline;
import fr.cnes.sirius.patrius.events.sensor.SensorVisibilityDetector;
import fr.cnes.sirius.patrius.frames.FramesFactory;
import fr.cnes.sirius.patrius.frames.TopocentricFrame;
import fr.cnes.sirius.patrius.math.util.FastMath;
import fr.cnes.sirius.patrius.orbits.Orbit;
import fr.cnes.sirius.patrius.propagation.BoundedPropagator;
import fr.cnes.sirius.patrius.propagation.analytical.KeplerianPropagator;
import fr.cnes.sirius.patrius.propagation.events.ConstantRadiusProvider;
import fr.cnes.sirius.patrius.propagation.events.EventDetector;
import fr.cnes.sirius.patrius.propagation.events.ThreeBodiesAngleDetector;
import fr.cnes.sirius.patrius.time.AbsoluteDate;
import fr.cnes.sirius.patrius.time.AbsoluteDateInterval;
import fr.cnes.sirius.patrius.utils.exception.PatriusException;
import fr.cnes.sirius.patrius.utils.exception.PropagationException;
import reader.Site;
import utils.ConstantsBE;
import utils.ProjectUtilities;
import utils.VTSTools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;


/**
 * This class encapsulates the context of an Earth Observation mission.
 *
 * @author herberl
 */
public class CompleteMission extends SimpleMission {

    /**
     * Maximum checking interval (s) for the event detection during the orbit
     * propagation.
     */
    public static final double MAXCHECK_EVENTS = 120.0;

    /**
     * Default convergence threshold (s) for the event computation during the orbit
     * propagation.
     */
    public static final double TRESHOLD_EVENTS = 1.e-4;

    /**
     * Arbitrary margin to ensure feasible slew
     */
    public static final double SLEW_MARGIN = 0.01;

    /**
     * This {@link Map} will be used to enumerate each site access {@link Timeline},
     * that is to say a {@link Timeline} with access windows respecting all
     * observation conditions. This object corresponds to the access plan, which
     * will be computed in the computeAccessPlan() method.
     */
    private final Map<Site, Timeline> accessPlan;

    /**
     * This {@link Map} will be used to enumerate each site's programmed
     * observation. We suggest to use an {@link AttitudeLawLeg} to encapsulate the
     * guidance law of each observation. This object corresponds to the observation
     * plan, which will be computed in the computeObservationPlan() method.
     */
    private final Map<Site, AttitudeLawLeg> observationPlan;

    /**
     * {@link StrictAttitudeLegsSequence} representing the cinematic plan during the
     * whole mission horizon. Each {@link AttitudeLeg} corresponds to a different
     * attitude law : either nadir pointing, target pointing or a slew between two
     * laws. This object corresponds to the cinematic plan, which will be computed
     * in the computeCinematicPlan() method.
     */
    private final StrictAttitudeLegsSequence<AttitudeLeg> cinematicPlan;

    /**
     * Constructor for the {@link CompleteMission} class.
     *
     * @param missionName   Name of the mission
     * @param numberOfSites Number of target {@link Site} to consider, please give a
     *                      number between 1 and 100.
     * @throws PatriusException      Can be raised by Patrius when building
     *                               particular objects. Here it comes from
     *                               {@link FramesFactory}
     * @throws IllegalStateException if the mission horizon is too short
     */
    public CompleteMission(final String missionName, int numberOfSites) throws PatriusException {

        // Since this class extends the SimpleMission class, we need to use the super
        // constructor to instantiate our instance of CompleteMission. All the
        // attributes of the super class will be instantiated during this step.
        super(missionName, numberOfSites);

        // Initialize the mission plans with empty maps. You will fill those HashMaps in
        // the "compute****Plan()" methods.
        this.accessPlan = new HashMap<>();
        this.observationPlan = new LinkedHashMap<>();
        this.cinematicPlan = new StrictAttitudeLegsSequence<>();

    }

    /**
     * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
     * <p>
     * Compute the access plan.
     * <p>
     * Reminder : the access plan corresponds to the object gathering all the
     * opportunities of access for all the sites of interest during the mission
     * horizon. One opportunity of access is defined by an access window (an
     * interval of time during which the satellite can observe the target and during
     * which all the observation conditions are achieved : visibility, incidence
     * angle, illumination of the scene,etc.). Here, we suggest you use the Patrius
     * class {@link Timeline} to encapsulate all the access windows of each site of
     * interest. One access window will then be described by the {@link Phenomenon}
     * object, itself defined by two {@link CodedEvent} objects giving the start and
     * end of the access window. Please find more tips and help in the submethods of
     * this method.
     *
     * @return the sites access plan with one {@link Timeline} per {@link Site}
     * @throws PatriusException If a {@link PatriusException} occurs during the
     *                          computations
     */
    public Map<Site, Timeline> computeAccessPlan() throws PatriusException {
        /**
         * Here you need to compute one access Timeline per target Site. You can start
         * with only one site and then try to compute all of them.
         * <p>
         * Note : when computing all the sites, try to make sure you don't decrease the
         * performance of the code too much. You might have some modifications to do in
         * order to ensure a reasonable time of execution.
         */
        /*
         * We give a very basic example of incomplete code computing the first target
         * site access Timeline and adding it to the accessPlan.
         *
         * Please complete the code below.
         */

        for (Site targetSite : this.getSiteList()) {
            Timeline siteAccessTimeline = createSiteAccessTimeline(targetSite);
            accessPlan.put(targetSite, siteAccessTimeline);
            ProjectUtilities.printTimeline(siteAccessTimeline);
        }
        /*
        Site targetSite = this.getSiteList().get(0);
        Timeline siteAccessTimeline = createSiteAccessTimeline(targetSite);
        accessPlan.put(targetSite, siteAccessTimeline);
        ProjectUtilities.printTimeline(siteAccessTimeline);
        */
        return this.accessPlan;
    }

    public Map<Site, AttitudeLawLeg> computeObservationPlan() throws PatriusException {
        
        // Clear the observation plan
        this.observationPlan.clear();

        // Initialization, nadir law
        AbsoluteDate previousObsEnd = ConstantsBE.START_DATE;
        AttitudeLaw previousObservationLaw = this.getSatellite().getDefaultAttitudeLaw();

        // Cities to observe by cloning the list of cities
        ArrayList<Site> remainingSites = new ArrayList<>(this.getSiteList());

        // Add observations while possible
        boolean cityToAdd = true;
        while (cityToAdd) {
        
            // Initialization of the best next city
            Site bestNextSite = null;
            double bestScore = Double.POSITIVE_INFINITY;
            AbsoluteDate bestObsStart = null;
            AbsoluteDate bestObsEnd = null;
            AttitudeLaw bestObservationLaw = null;

            // Loop on all the remaining cities to find the best next city (minimize slew_time)
            for (Site nextTarget : remainingSites) {

                // Find the observation beginning for this city candidate
                AttitudeLaw nextTargetObservationLaw = this.createObservationLaw(nextTarget);

                // Getting its access Timeline
                final Timeline timeline = this.accessPlan.get(nextTarget);
                AbsoluteDate obsStart = null;
                AbsoluteDate obsEnd = null;

                // Loop on all the access windows to find the first one that is valid
                for (final Phenomenon accessWindow : timeline.getPhenomenaList()) {
                    // The Phenomena are sorted chronologically so the accessIntervals List is too
                    AbsoluteDateInterval accessInterval = accessWindow.getTimespan();

                    // Getting the beginning/end of the accessInterval as AbsoluteDate objects
                    AbsoluteDate accessStart = accessInterval.getLowerData();
                    AbsoluteDate accessEnd = accessInterval.getUpperData();

                    if (accessEnd.durationFrom(accessStart) < ConstantsBE.INTEGRATION_TIME) {
                        // If the access window is too short, we skip it
                        continue;
                    }

                    if (accessEnd.compareTo(previousObsEnd) < 0) {
                        // The access window ends before the previous observation ends, we skip it
                        continue;
                    }

                    obsStart = previousObsEnd;
                    if (accessStart.durationFrom(previousObsEnd) > this.getSatellite().getMaxSlewDuration()) {
                        // If the next access window starts far enough from the previous observation to perform any slew, we use
                        // the beginning of the next window  for the observation.
                        obsStart = accessStart;
                    } else { // Otherwise, we process the start of the observation with the slew duration
                        ConstantSpinSlew slew = this.computeExactSlew(previousObservationLaw, nextTargetObservationLaw,
                                previousObsEnd, true, "");
                        double slewDuration = slew.getDuration();

                        if (accessEnd.durationFrom(previousObsEnd) < ConstantsBE.INTEGRATION_TIME + slewDuration) {
                            // If the access window ends before we can reach it by slewing and observe it during the full duration, we skip it
                            continue;
                        }

                        // Observation start is the end of the slew
                        obsStart = slew.getEnd();
                        // If access window starts after the end of the slew, we use the beginning of the access window
                        if (accessStart.compareTo(obsStart) > 0) {
                            obsStart = accessStart;
                        }
                    }

                    // Process the end of the observation and stop the loop because earliest observation found
                    obsEnd = obsStart.shiftedBy(ConstantsBE.INTEGRATION_TIME);
                    break;
                }
                // At this point, we found the earliest possible observation start for this city candidate
                // Now we compare that candidate to the best one found so far
                if (obsStart == null || obsEnd == null) {
                    // If the observation start/end are null, it means that we could not find a valid observation for this city candidate
                    continue;
                }
                // Compute the score of the candidate and update the best candidate if needed
                double nextTargetScore = obsStart.durationFrom(previousObsEnd);
                if (nextTargetScore < bestScore) {
                    bestNextSite = nextTarget;
                    bestScore = nextTargetScore;
                    bestObsStart = obsStart;
                    bestObsEnd = obsEnd;
                    bestObservationLaw = nextTargetObservationLaw;
                }
            }

            // We found the best observation to add to the plan
            // Check if it's not null (no more city to add to the plan)
            if (bestNextSite == null) {
                cityToAdd = false;
                continue;
            }

            // Else we add it to the plan
            String legName = "OBS_" + bestNextSite.getName();
            AbsoluteDateInterval obsInterval = new AbsoluteDateInterval(bestObsStart, bestObsEnd);
            AttitudeLawLeg obsLeg = new AttitudeLawLeg(bestObservationLaw, obsInterval, legName);
            this.observationPlan.put(bestNextSite, obsLeg);

            // The city is added, remove it from the list of remaining cities
            remainingSites.remove(bestNextSite);

            // Print the observation plan
            System.out.println("Selected observation : " + bestNextSite.getName() + " from " + bestObsStart + " to " + bestObsEnd);

            // And we update the variables for the next iteration
            previousObsEnd = bestObsEnd;
            previousObservationLaw = bestObservationLaw;
        }
        return this.observationPlan;
    }

    /**
     * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
     * <p>
     * Computes the cinematic plan.
     * <p>
     * Here you need to compute the cinematic plan, which is the cinematic chain of
     * attitude law legs (observation, default law and slews) needed to perform the
     * mission. Usually, we start and end the mission in default law and during the
     * horizon, we alternate between default law, observation legs and slew legs.
     *
     * @return a {@link StrictAttitudeLegsSequence} that gives all the cinematic
     * plan of the {@link Satellite}. It is a chronological sequence of all
     * the {@link AttitudeLawLeg} that are necessary to define the
     * {@link Attitude} of the {@link Satellite} during all the mission
     * horizon. Those legs can have 3 natures : pointing a target site,
     * pointing nadir and performing a slew between one of the two previous
     * kind of legs.
     * @throws PatriusException in case of issue
     */
    public StrictAttitudeLegsSequence<AttitudeLeg> computeCinematicPlan() throws PatriusException {

        /**
         * Now we want to assemble a continuous attitude law which is valid during all
         * the mission horizon. For that, we will use to object
         * StrictAttitudeLegsSequence<AttitudeLeg> which is a chronological sequence of
         * AttitudeLeg. In our case, each AttitudeLeg will be an AttitudeLawLeg, either
         * a leg of site observation, a slew, or the nadir pointing attitude law (see
         * the Satellite constructor and the BodyCenterGroundPointing class, it's the
         * Satellite default attitude law). For more help about the Attitude handling,
         * use the module 11 of the patrius formation.
         * <p>
         * Tip 1 : Please give names to the different AttitudeLawLeg you build so that
         * you can visualize them with VTS later on. For example "OBS_Paris" when
         * observing Paris or "SlEW_Paris_Lyon" when adding a slew from Paris
         * observation AttitudeLawLeg to Lyon observation AttitudeLawLeg.
         * <p>
         * Tip 2 : the sequence you want to obtain should look like this :
         * [nadir-slew-obs1-slew-obs2-slew-obs3-slew-nadir] for the simple version where
         * you don't try to fit nadir laws between observations or
         * [nadir-slew-obs1-slew-nadir-slew-obs2-slew-obs3-slew-nadir] for the more
         * complex version with nadir laws if the slew during two observation is long
         * enough.
         * <p>
         * Tip 3 : You can use the class ConstantSpinSlew(initialAttitude,
         * finalAttitude, slewName) for the slews. This an AttitudeLeg so you will be
         * able to add it to the StrictAttitudeLegsSequence as every other leg.
         */

        /*
         * Example of code using our observation plan, let's say we only have one obs
         * pointing Paris.
         *
         * Then we are going to create a very basic cinematic plan : nadir law => slew
         * => obsParis => slew => nadir law
         *
         * To do that, we need to compute the slew duration from the end of nadir law to
         * the beginning of Paris obs and then from the end of Paris obs to the beginning
         * of nadir law. For that, we use the Satellite#computeSlewDurationMethod() as
         * before. We know we have to the time to perform the slew thanks to the
         * cinematic checks we already did during the observation plan computation.
         */

        // Clear the cinematic plan
        this.cinematicPlan.clear();

        // Getting our nadir law
        final AttitudeLaw nadirLaw = this.getSatellite().getDefaultAttitudeLaw();
        AttitudeLeg previousAttitudeLeg = null;
        AbsoluteDate previousAttitudeDate = ConstantsBE.START_DATE;
        String previousTargetName = "";
        boolean isBeginning = true;


        // Scrolling through the entries of the accessPlan
        for (final Entry<Site, AttitudeLawLeg> entry : this.observationPlan.entrySet()) {
            // Getting the target Site
            final Site target = entry.getKey();

            // Getting the observation law of the target
            AttitudeLeg targetAttitudeLawLeg = this.observationPlan.get(target);
            AbsoluteDate obsStart = targetAttitudeLawLeg.getDate();
            AbsoluteDate obsEnd = targetAttitudeLawLeg.getEnd();


            // If at the beginning of the propagation, we are in nadir law
            if (isBeginning) {
                // Compute the slew between the nadir attitude at the start of the observation
                // and the target pointing attitude at the start of the slew
                ConstantSpinSlew slewToTarget = this.computeExactSlew(nadirLaw, targetAttitudeLawLeg,
                        obsStart, false, "NADIR_TO_" + target.getName());

                // Create the nadir attitude leg from the start of the propagation to the beginning of the first slew
                AbsoluteDateInterval nadirInterval = new AbsoluteDateInterval(
                        ConstantsBE.START_DATE, slewToTarget.getDate());
                AttitudeLawLeg nadirLeg = new AttitudeLawLeg(nadirLaw, nadirInterval, "START_NADIR");

                // Add both legs to the cinematic plan
                this.cinematicPlan.add(nadirLeg);
                this.cinematicPlan.add(slewToTarget);
                isBeginning = false;
            }

            else if (obsStart.durationFrom(previousAttitudeDate) > 3 * this.getSatellite().getMaxSlewDuration()) {
                // If the interval since the previous observation is long enough, we add a nadir law between them

                // Create the slew from the previous target to nadir
                ConstantSpinSlew slewToNadir = this.computeExactSlew(previousAttitudeLeg, nadirLaw,
                        previousAttitudeDate, true, previousTargetName + "_TO_NADIR");

                // Create the slew from nadir to the next target
                ConstantSpinSlew slewFromNadir = this.computeExactSlew(nadirLaw, targetAttitudeLawLeg,
                        obsStart, false, "NADIR_TO_" + target.getName());

                // Create the nadir attitude leg from the end of the first slew to the start of the second slew
                AbsoluteDateInterval nadirInterval = new AbsoluteDateInterval(
                        slewToNadir.getEnd(), slewFromNadir.getDate());
                AttitudeLawLeg nadirLeg = new AttitudeLawLeg(nadirLaw, nadirInterval,
                        "NADIR_BETWEEN_" + previousTargetName + "_AND_" + target.getName());

                // Add all three legs to the cinematic plan
                this.cinematicPlan.add(slewToNadir);
                this.cinematicPlan.add(nadirLeg);
                this.cinematicPlan.add(slewFromNadir);
            }

            else {
                // Slewing to the target
                Attitude attitude1 = previousAttitudeLeg.getAttitude(this.getSatellite().getPropagator().getPvProvider(), previousAttitudeDate, this.getEme2000());
                Attitude attitude2 = targetAttitudeLawLeg.getAttitude(this.getSatellite().getPropagator().getPvProvider(), obsStart, this.getEme2000());
                ConstantSpinSlew slewToTarget = new ConstantSpinSlew(attitude1, attitude2, previousTargetName + "_TO_" + target.getName());
                this.cinematicPlan.add(slewToTarget);
            }

            // Add the observation law to the cinematic plan
            this.cinematicPlan.add(targetAttitudeLawLeg);

            // Update the previous attitude leg and date
            previousAttitudeLeg = targetAttitudeLawLeg;
            previousAttitudeDate = obsEnd;
            previousTargetName = target.getName();
        }

        // Adding the last leg to nadir if needed
        if (previousAttitudeLeg != null) {
            ConstantSpinSlew slewToNadir = this.computeExactSlew(previousAttitudeLeg, nadirLaw,
                    previousAttitudeDate, true, previousTargetName + "_TO_NADIR");
            this.cinematicPlan.add(slewToNadir);
            previousAttitudeDate = slewToNadir.getEnd();
        }
        
        // Create the nadir attitude leg from the end of the last slew to end of the propagation
        AbsoluteDateInterval nadirInterval = new AbsoluteDateInterval(previousAttitudeDate, ConstantsBE.END_DATE);
        AttitudeLawLeg nadirLeg = new AttitudeLawLeg(nadirLaw, nadirInterval, "END_NADIR");
        this.cinematicPlan.add(nadirLeg);

        return this.cinematicPlan;
    }

    /**
     * [DO NOT MODIFY THIS METHOD]
     * <p>
     * Checks the cinematic plan and prints if it's ok or not.
     * <p>
     * We provide this method so that you can check if your cinematic plan doesn't
     * violate a cinematic constraint. It returns a boolean saying if the plan is
     * valid or not.
     *
     * @return A boolean indicating the cinematic validity of the plan.
     * @throws PatriusException if an error occurs during propagation
     */
    public boolean checkCinematicPlan() throws PatriusException {

        final KeplerianPropagator propagator = createDefaultPropagator();

        // Checking the cinematic plan validity
        boolean valid = true;
        for (final AttitudeLeg slew : this.cinematicPlan) {
            System.out.println(slew.getNature());

            final Attitude endAtt = slew.getAttitude(propagator, slew.getEnd(), this.getEme2000());
            final Attitude startAtt = slew.getAttitude(propagator, slew.getDate(), this.getEme2000());

            final boolean condition = slew.getTimeInterval().getDuration() > this.getSatellite()
                    .computeSlewDuration(startAtt, endAtt);

            if (condition) {
                System.out.println("Cinematic is ok");
            } else {
                valid = false;
                System.out.println("WARNING : cinematic is not realist for this slew");
                System.out.println("Slew actual duration : " + slew.getTimeInterval().getDuration());
                System.out
                        .println("Slew duration theory : " + this.getSatellite().computeSlewDuration(startAtt, endAtt));
            }
        }
        System.out.println("==== Is the cinematic plan valid ? => " + valid + " ==== ");
        return valid;
    }

    /**
     * [DO NOT MODIFY THIS METHOD]
     * <p>
     * Compute the final score from the observation plan.
     *
     * <p>
     * Note : the observation plan should have unique sites.<br>
     * The duplicated elements won't be considered for the score (target with no
     * value = lost opportunity)
     * </p>
     *
     * @return the final score
     */
    public double computeFinalScore() {

        // Convert the observation plan into a Set of sites to make sure a site will not
        // be evaluated twice
        // The Set makes sure there won't be duplicated sites
        final Set<Site> sitesSet = new HashSet<>(this.observationPlan.size());
        for (final Entry<Site, AttitudeLawLeg> entry : this.observationPlan.entrySet()) {
            sitesSet.add(entry.getKey());
        }

        // Loop over each site and sum its score
        double finalScore = 0.;
        for (final Site site : sitesSet) {
            finalScore += site.getScore();
        }

        return finalScore;
    }

    /**
     * [DO NOT MODIFY THIS METHOD]
     * <p>
     * Writes the VTS output files : one CIC-POI file to print the sites of
     * interest, one CIC-OEM file giving the position and velocity ephemeris of the
     * satellite, one CIC-AEM file giving the attitude ephemeris of the satellite
     * pointing Nadir only (to help visualize the access field of view of the
     * satellite) and one CIC-AEM file giving the attitude ephemeris of the
     * satellite cinematic plan. Also writes the cinematic plan as a sequence of
     * pointing modes for the satellite in a CIC-MEM file.
     *
     * @param cinematicPlan Input cinematic plan.
     * @throws PropagationException if an error happens during the {@link Orbit}
     *                              propagation
     */
    public void generateVTSVisualization(StrictAttitudeLegsSequence<AttitudeLeg> cinematicPlan)
            throws PropagationException {
        // First, create the propagator for the satellite's pointing capacity view
        // (nadir law)
        final KeplerianPropagator vtsPropagatorNadir = createDefaultPropagator();
        vtsPropagatorNadir.setEphemerisMode();
        vtsPropagatorNadir.propagate(this.getEndDate());

        // Get generated ephemeris
        final BoundedPropagator ephemerisNadir = vtsPropagatorNadir.getGeneratedEphemeris();

        // Then, we create the propagator for the cinematic plan visualization
        final KeplerianPropagator vtsPropagator = new KeplerianPropagator(this.getSatellite().getInitialOrbit(),
                cinematicPlan);

        vtsPropagator.setEphemerisMode();
        vtsPropagator.propagate(this.getEndDate());

        // Get generated ephemeris
        final BoundedPropagator ephemeris = vtsPropagator.getGeneratedEphemeris();

        // Writing the outputs
        final String pathPOI = ConstantsBE.PATH_VTS_DIRECTORY + File.separator + "BE_Supaero_Target_Sites_POI.txt";
        final String pathOEM = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
                + "BE_Supaero_Satellite_Trajectory_OEM.txt";
        final String pathAEMNadir = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
                + "BE_Supaero_Nadir_Pointing_AEM.txt";
        final String pathAEMCinematicPlan = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
                + "BE_Supaero_Cinematic_Plan_AEM.txt";
        final String pathMEMCinematicPlan = ConstantsBE.PATH_VTS_DIRECTORY + File.separator
                + "BE_Supaero_Cinematic_Plan_Events_MEM.txt";

        System.out.println("\n\nWriting VTS outputs, please wait...");
        VTSTools.generatePOIFile(pathPOI, this.getSiteList());
        VTSTools.generateOEMFile(pathOEM, this.getStartDate(), this.getEndDate(), ephemeris);
        VTSTools.generateAEMFile(pathAEMNadir, this.getStartDate(), this.getEndDate(), ephemerisNadir);
        VTSTools.generateAEMFile(pathAEMCinematicPlan, this.getStartDate(), this.getEndDate(), ephemeris);
        VTSTools.generateLegSequenceMEMFile(pathMEMCinematicPlan, cinematicPlan);
        System.out.println("VTS outputs written");
    }

    /**
     * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
     * <p>
     * This method should compute the input {@link Site}'s access {@link Timeline}.
     * That is to say the {@link Timeline} which contains all the {@link Phenomenon}
     * respecting the access conditions for this site : good visibility + correct
     * illumination of the {@link Site}.
     * <p>
     * For that, we suggest you create as many {@link Timeline} as you need and
     * combine them with logical gates to filter only the access windows phenomenon.
     *
     * @param targetSite Input target {@link Site}
     * @return The {@link Timeline} of all the access {@link Phenomenon} for the
     * input {@link Site}.
     * @throws PatriusException If a {@link PatriusException} occurs.
     */
    private Timeline createSiteAccessTimeline(Site targetSite) throws PatriusException {
        /**
         * Step 1 :
         * <p>
         * Create one Timeline per phenomenon you want to monitor.
         */
        /*
         * Use the createSiteXTimeline method to create a custom Timeline. More
         * indication are given inside the method. Note that you will have to code one
         * method per constraint, for example the method createSiteVisibilityTimeline
         * for visibility constraint and createSiteIlluminationTimeline for illumination
         * constraint. All the methods you code can be coded using the given
         * createSiteXTimeline method as a basis.
         */
        Timeline timeline1 = createSiteVisibilityTimeline(targetSite);
        Timeline timeline2 = createSiteIlluminationTimeline(targetSite);
        Timeline timeline3 = createSitePhaseAngleTimeline(targetSite);
        // etc.

        /**
         * Step 2 :
         * <p>
         * Combine the timelines with logical gates and retrieve only the access
         * conditions through a refined Timeline object.
         * <p>
         * For that, you can use the classes in the events.postprocessing module : for
         * example, the AndCriterion or the NotCriterion.
         * <p>
         * Finally, you can filter only the Phenomenon matching a certain condition
         * using the ElementTypeFilter
         */
        /*
         * Code your logical operations on Timeline objects and filter only the access
         * Phenomenon (gathering all constraints you need to define an access condition)
         * below.
         */
        // Combining all Timelines
        // Creating a global Timeline containing all phenomena, this Timeline will be
        // filtered and processed to that only the access Phenomenon remain, this is
        // our siteAccessTimeline
        final Timeline siteAccessTimeline = new Timeline(
                new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()));
        // Adding the phenomena of all the considered timelines
        for (final Phenomenon phenom : timeline1.getPhenomenaList()) {
            siteAccessTimeline.addPhenomenon(phenom);
        }
        for (final Phenomenon phenom : timeline2.getPhenomenaList()) {
            siteAccessTimeline.addPhenomenon(phenom);
        }
        for (final Phenomenon phenom : timeline3.getPhenomenaList()) {
            siteAccessTimeline.addPhenomenon(phenom);
        }

        // Define and use your own criteria, here is an example (use the right strings
        // defined when naming the phenomenon in the GenericCodingEventDetector)
        // Applying our criterion adds all the new phenomena inside the global timeline
        AndCriterion andCriterion = new AndCriterion("Visibility", "Illumination",
                "Visibility AND Illumination", "OBS OK");
        andCriterion.applyTo(siteAccessTimeline);

        AndCriterion andCriterion2 = new AndCriterion("Visibility AND Illumination", "Phase angle",
                "Visibility AND Illumination AND Phase angle", "OBS OK (BONUS)");
        andCriterion2.applyTo(siteAccessTimeline);


        // Then create an ElementTypeFilter that will filter all phenomenon not
        // respecting the input condition you gave it
        final ElementTypeFilter obsConditionFilter = new ElementTypeFilter("Visibility AND Illumination AND Phase angle", false);
        // Finally, we filter the global timeline to keep only X1 AND X2 phenomena
        obsConditionFilter.applyTo(siteAccessTimeline);

        /*
         * Now make sure your globalTimeline represents the access Timeline for the
         * input target Site, and it's done ! You can print the Timeline using the
         * utility module of the BE as below
         */

        // Log the final access timeline associated to the current target
        System.out.println("\n" + targetSite.getName());
        // ProjectUtilities.printTimeline(siteAccessTimeline);

        return siteAccessTimeline;
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return The {@link Timeline} containing all the {@link Phenomenon} relative
     * to the visibility of the target.
     * @throws PatriusException If a {@link PatriusException} occurs when creating
     *                          the {@link Timeline}.
     */
    private Timeline createSiteVisibilityTimeline(Site targetSite) throws PatriusException {
        // Building the EventDetector
        EventDetector constraintVisibilityDetector = createConstraintVisibilityDetector(targetSite);

        // Adding the detector to the Orbit Propagator
        KeplerianPropagator propagator = this.createDefaultPropagator();
        propagator.addEventDetector(constraintVisibilityDetector);

        // Creating the GenericCodingEventDetector and the CodedEventsLogger to detect the events and visualize them
        GenericCodingEventDetector codingEventVisibilityDetector = new GenericCodingEventDetector(constraintVisibilityDetector,
                "Entering visibility", "Leaving visibility", true, "Visibility");
        CodedEventsLogger eventVisibilityLogger = new CodedEventsLogger();
        EventDetector eventVisibilityDetector = eventVisibilityLogger.monitorDetector(codingEventVisibilityDetector);
        // Adding the logger to the propagator to monitor the event coded by the codingEventDetector
        propagator.addEventDetector(eventVisibilityDetector);


        // Propagating the orbit
        propagator.propagate(this.getStartDate(), this.getEndDate());

        // Post-processing the events: creating a Timeline containing all the visibility Phenomenon
        return new Timeline(eventVisibilityLogger,
                new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), this.getSatellite().getSpacecraftState());
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return The {@link Timeline} containing all the {@link Phenomenon} relative
     * to the illumination of the target.
     * @throws PatriusException If a {@link PatriusException} occurs when creating
     *                          the {@link Timeline}.
     */
    private Timeline createSiteIlluminationTimeline(Site targetSite) throws PatriusException {

        // Building the EventDetector
        EventDetector constraintIlluminationDetector = createConstraintIlluminationDetector(targetSite);

        // Adding the detector to the Orbit Propagator
        KeplerianPropagator propagator = this.createDefaultPropagator();
        propagator.addEventDetector(constraintIlluminationDetector);

        // Creating the GenericCodingEventDetector and the CodedEventsLogger to detect the events and visualize them
        GenericCodingEventDetector codingEventIlluminationDetector =
                new GenericCodingEventDetector(constraintIlluminationDetector, "Leaving illumination",
                        "Entering illumination", false, "Illumination");
        CodedEventsLogger eventIlluminationLogger = new CodedEventsLogger();
        EventDetector eventIlluminationDetector = eventIlluminationLogger.monitorDetector(codingEventIlluminationDetector);
        // Adding the logger to the propagator to monitor the event coded by the codingEventDetector
        propagator.addEventDetector(eventIlluminationDetector);

        // Propagating the orbit
        propagator.propagate(this.getStartDate(), this.getEndDate());

        // Post-processing the events: creating a Timeline containing all the illumination Phenomenon
        return new Timeline(eventIlluminationLogger,
                new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), this.getSatellite().getSpacecraftState());
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return The {@link Timeline} containing all the {@link Phenomenon} relative
     * to the blinding of the sensor, which is defined by the phase angle relative to the Sun.
     * @throws PatriusException If a {@link PatriusException} occurs when creating
     *                          the {@link Timeline}.
     */
    private Timeline createSitePhaseAngleTimeline(Site targetSite) throws PatriusException {

        // Building the EventDetector
        EventDetector constraintPhaseAngleDetector = createConstraintPhaseAngleDetector(targetSite);

        // Adding the detector to the Orbit Propagator
        KeplerianPropagator propagator = this.createDefaultPropagator();
        propagator.addEventDetector(constraintPhaseAngleDetector);


        // Creating the GenericCodingEventDetector and the CodedEventsLogger to detect the events and visualize them
        GenericCodingEventDetector codingEventPhaseAngleDetector =
                new GenericCodingEventDetector(constraintPhaseAngleDetector, "Sun blinding sensor",
                        "Sun behind sensor", false, "Phase angle");
        CodedEventsLogger eventPhaseAngleLogger = new CodedEventsLogger();
        EventDetector eventPhaseAngleDetector = eventPhaseAngleLogger.monitorDetector(codingEventPhaseAngleDetector);

        // Adding the logger to the propagator to monitor the event coded by the codingEventDetector
        propagator.addEventDetector(eventPhaseAngleDetector);

        // Propagating the orbit
        propagator.propagate(this.getStartDate(), this.getEndDate());

        // Post-processing the events: creating a Timeline containing all the sensor blinding Phenomenon
        return new Timeline(eventPhaseAngleLogger,
                new AbsoluteDateInterval(this.getStartDate(), this.getEndDate()), this.getSatellite().getSpacecraftState());
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return An {@link EventDetector} answering the visibility constraint
     */
    private EventDetector createConstraintVisibilityDetector(Site targetSite) {
        SensorModel sensorModel = new SensorModel(this.getSatellite().getAssembly(), Satellite.SENSOR_NAME);
        sensorModel.addMaskingCelestialBody(this.getEarth());
        sensorModel.setMainTarget(new TopocentricFrame(this.getEarth(), targetSite.getPoint(), targetSite.getName()),
                new ConstantRadiusProvider(0.0)); // Radius = 0 --> assumption of point target

        return new SensorVisibilityDetector(sensorModel, MAXCHECK_EVENTS, TRESHOLD_EVENTS,
                EventDetector.Action.CONTINUE, EventDetector.Action.CONTINUE);
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return An {@link EventDetector} answering the illumination constraint
     */
    private EventDetector createConstraintIlluminationDetector(Site targetSite) {
        return new ThreeBodiesAngleDetector(new TopocentricFrame(this.getEarth(), targetSite.getPoint(),
                targetSite.getName()), this.getEarth(), this.getSun(),
                FastMath.toRadians(ConstantsBE.MAX_SUN_INCIDENCE_ANGLE),
                MAXCHECK_EVENTS, TRESHOLD_EVENTS, EventDetector.Action.CONTINUE);
    }

    /**
     * @param targetSite Input target {@link Site}
     * @return An {@link EventDetector} answering the sensor blinding constraint
     */
    private EventDetector createConstraintPhaseAngleDetector(Site targetSite) {
        return new ThreeBodiesAngleDetector(this.getSatellite().getPropagator().getPvProvider(), new TopocentricFrame(this.getEarth(), targetSite.getPoint(),
                targetSite.getName()), this.getSun(), FastMath.toRadians(ConstantsBE.MAX_SUN_PHASE_ANGLE),
                MAXCHECK_EVENTS, TRESHOLD_EVENTS, EventDetector.Action.CONTINUE);
    }

    /**
     * [COMPLETE THIS METHOD TO ACHIEVE YOUR PROJECT]
     * <p>
     * Create an observation leg, that is to say an {@link AttitudeLaw} that give
     * the {@link Attitude} (pointing direction) of the {@link Satellite} in order
     * to perform the observation of the input target {@link Site}.
     * <p>
     * An {@link AttitudeLaw} is an {@link AttitudeProvider} providing the method
     * {@link AttitudeProvider#getAttitude(Orbit)} which can be used to compute the
     * {@link Attitude} of the {@link Satellite} at any given {@link AbsoluteDate}
     * (instant) during the mission horizon.
     * <p>
     * An {@link AttitudeLaw} is valid at anu time in theory.
     *
     * @param target Input target {@link Site}
     * @return An {@link AttitudeLawLeg} adapted to the observation.
     */
    private AttitudeLaw createObservationLaw(Site target) {
        /**
         * To perform an observation, the satellite needs to point the target for a
         * fixed duration.
         * <p>
         * Here, you will use the {@link TargetGroundPointing}. This law provides the
         * Attitude of a Satellite that only points one target at the surface of a
         * BodyShape. The earth object from the SimpleMission is a BodyShape and we
         * remind you that the Site object has an attribute which is a GeodeticPoint.
         * Use that information to your advantage to build a TargetGroundPointing.
         */
        /*
         * Complete the code below to create your observation law and return it
         */
        return new TargetGroundPointing(
            this.getEarth(), target.getPoint(), this.getSatellite().getSensorAxis(),
            this.getSatellite().getFrameXAxis());
    }

    /**
     * Process the slew between two attitude laws while the satellite is moving.
     *
     * @param attitudeProvider1 First attitude law
     * @param attitudeProvider2 Second attitude law
     * @param date              Fixed date constraint (start or end of the slew)
     * @param forward           Boolean indicating if the date is the start or the end of the slew
     * @param name              Name of the slew
     * @return A {@link Timeline} containing the slew
     * @throws PatriusException If a {@link PatriusException} occurs when creating
     */
    private ConstantSpinSlew computeExactSlew(AttitudeProvider attitudeProvider1, AttitudeProvider attitudeProvider2,
                                              AbsoluteDate date, boolean forward, String name) throws PatriusException {
        // initialization of the variables
        double error = Double.POSITIVE_INFINITY;
        double previousSlewDuration = 0.0;
        double slewDuration = 0.0;
        AbsoluteDate beginDate = date;
        AbsoluteDate endDate = date;

        // process the attitudes at the beginning and the end of the slew
        Attitude attitude1 = attitudeProvider1.getAttitude(
                this.getSatellite().getPropagator().getPvProvider(), beginDate, this.getEme2000());
        Attitude attitude2 = attitudeProvider2.getAttitude(
                this.getSatellite().getPropagator().getPvProvider(), endDate, this.getEme2000());
        
        // while error in slew duration is greater than tolerance
        while (error > CompleteMission.SLEW_MARGIN) {
                // compute the slew duration between the two attitudes
                slewDuration = this.getSatellite().computeSlewDuration(attitude1, attitude2);

                // if forward, update the end date, else update the beginning date
                if (forward) {
                    endDate = beginDate.shiftedBy(slewDuration);
                } else {
                    beginDate = endDate.shiftedBy(-slewDuration);
                }

                // update the error and the previous slew duration
                error = 10*Math.abs(slewDuration - previousSlewDuration);
                previousSlewDuration = slewDuration;

                // update the attitudes
                attitude1 = attitudeProvider1.getAttitude(this.getSatellite().getPropagator().getPvProvider(), beginDate, this.getEme2000());
                attitude2 = attitudeProvider2.getAttitude(this.getSatellite().getPropagator().getPvProvider(), endDate, this.getEme2000());
        }
        // add tolerance to the adaptable date to ensure feasibility
        if (forward) {
                endDate = endDate.shiftedBy(CompleteMission.SLEW_MARGIN);
        } else {
                beginDate = beginDate.shiftedBy(-CompleteMission.SLEW_MARGIN);
        }
        // process the final attitudes
        attitude1 = attitudeProvider1.getAttitude(this.getSatellite().getPropagator().getPvProvider(), beginDate, this.getEme2000());
        attitude2 = attitudeProvider2.getAttitude(this.getSatellite().getPropagator().getPvProvider(), endDate, this.getEme2000());
        // create the final slew
        return new ConstantSpinSlew(attitude1, attitude2, name);
    }

    @Override
    public String toString() {
        return "CompleteMission [name=" + this.getName() + ", startDate=" + this.getStartDate() + ", endDate="
                + this.getEndDate() + ", satellite=" + this.getSatellite() + "]";
    }
}
