# Teaching Earth science with Mantle-Bloom

Mantle-Bloom runs an Earth-like planet forward in time and lets learners inspect how its
terrain, climate, and water systems interact. A classroom activity can use those views to
make a process visible, test an explanation, or explore the consequences of a change.
The research below comes from other learning environments. The proposed pilot would test
how these ideas work with Mantle-Bloom. [Simulation model](../docs/simulation-model.md)

## When the model can help students learn

### Make the evidence easy to find

A mystery landscape gives students a reason to look for evidence. Give them a marked
region, two plausible explanations, and a worked example of reading one layer. Ask them
to connect two observations to a process, then revise a sketch or written explanation
after discussing the evidence. In a study of 950 secondary students using Tectonic
Explorer, 7% gave a fully sequenced causal explanation for a convergence task and 3%
did so for a harder task involving processes on opposite sides of a plate. Another
middle-school study found that students with stronger spatial skills were more likely
to interpret depth and notice patterns in a visualization. A grade-eight classroom
study offers a precedent for pairing discussion with student-generated drawings.
[Pallant et al., 2023](https://doi.org/10.1080/09500693.2023.2217471),
[Focus on the notice](https://link.springer.com/article/10.1186/s41235-020-00263-0),
[Grade-eight drawing study](https://www.tandfonline.com/doi/full/10.1080/09500693.2023.2214688)

### Use playback when the sequence or rate matters

For a deep-time lesson, pause before a change and ask students what will move next.
Compare dated frames using the same map view and a distance scale. Students can estimate
how quickly a landmark moves and compare their estimate with measured plate-motion
rates. For a question about where precipitation falls, labeled stills can give students
time to read the maps. A review of instructional animations found the clearest learning
benefit when students needed to understand the features of change itself.
[Ploetzner et al., 2020](https://onlinelibrary.wiley.com/doi/10.1111/jcal.12476),
[USGS Cascadia GPS data](https://www.usgs.gov/special-topics/subduction-zone-science/science/cascadia-subduction-zone-database)

### Compare one change and explain the outcome

Start a comparison from the same saved state, change one setting, and measure one
outcome at the same elapsed time. Ask for a prediction and an explanation of the observed
difference. A high-school extension could compare precipitation over the same terrain
using Mantle-Bloom's diagnostic and CFD wind settings. Students could identify one
conclusion supported by both runs and one that depends on the setting. The model
documentation reports substantial precipitation differences in one such comparison.
[Mantle-Bloom wind model](../docs/simulation-model/climate.md#wind-model-diagnostic-default-vs-cfd)

Start a systems lesson with a few causal links students can explain. For example, trace
how a range changes precipitation and how precipitation can affect erosion. A follow-up
can ask how the resulting change in terrain might alter precipitation again. In an
eighth-grade water-cycle study, many students gained a more connected view of the system.
Reasoning about hidden processes and cycles remained harder.
[Assaraf and Orion, 2005](https://onlinelibrary.wiley.com/doi/10.1002/tea.20061)

### Give teachers prompts, examples, and a way to check transfer

Teachers in the GEODE plate tectonics study particularly valued sample student answers
and explanations of why answers were correct or incorrect. A three-year study of two
high-school teachers documented how they adapted a model-based climate curriculum to
their students and schedules. For Mantle-Bloom, a one-page teacher guide could include
a marked map, a sample explanation, a common incomplete explanation, and a discussion
prompt. An extension could add a quantitative comparison and a question about model
assumptions. [Lord et al., teacher study](https://concord.org/wp-content/uploads/publications/a-remote-view-into-the-classroom.pdf),
[Carroll Steward et al., high-school climate curriculum](https://www.giss.nasa.gov/pubs/abs/ca05030t.html)

A meta-analysis of 72 inquiry-learning studies found a positive average effect of
guidance on learning outcomes. A useful lesson sequence is to ask students to notice
a relevant difference, explain it, get feedback, and apply the explanation to a fresh
case. The fresh case gives the teacher a way to assess whether students can use the
idea independently. [Lazonder and Harmsen, 2016](https://journals.sagepub.com/doi/10.3102/0034654315627366)

## A possible first pilot

A good first pilot might give students in grades 7–9 one class period to investigate
**why one side of a mountain range is drier**. The learning goal is to use evidence
about prevailing wind and terrain to explain a precipitation pattern.

1. Show a marked range on a simulated world and ask students to predict the wetter side
   from the terrain and wind. Ask each student to record a reason.
2. Reveal the model's precipitation map. Have pairs annotate where moist air rises and
   cools, where precipitation falls, and why less moisture reaches the far side.
3. Compare the simulated pattern with observed precipitation across Washington's
   Cascades. Ask students to point to the evidence for their explanation on both maps.
4. Give each student a fresh map with wind coming from the opposite direction. Score
   whether they identify the wetter side, explain how wind and terrain produce the
   pattern, and cite a feature of the map as evidence.

Prepared frames and a one-page teacher guide are enough to try the sequence. Record
which steps fit in the period and where students' explanations break down. A later
comparison could give two groups the same prompts and transfer question, with one
group using prepared frames and the other using the interactive model.
[Washington State Climate Office](https://climate.uw.edu/washingtons-climate/),
[PRISM precipitation normals](https://www.prism.oregonstate.edu/normals/)

## Other uses to explore

- **Mystery planet.** Students reconstruct a landscape's history from clues such as
  plate boundaries, relief, and river networks. A teacher-controlled reveal could
  expose earlier checkpoints after students propose competing explanations.
- **Sketch-driven experiments.** Students sketch a continent with a mountain barrier
  and compare it with a version whose barrier sits at a different latitude. To make
  a causal claim, they would need to check which other world properties changed.
  [Sketch-driven world generation](../docs/simulation-model/foundations.md#human-made-worlds-sketch-driven-generation-worldsketchpy)
- **Ice and coastlines.** Students compare dated states from a glacial cycle and use
  ice and shoreline evidence to explain changing sea level. An extension could compare
  grounded and floating ice using a water-budget diagram.
  [Ice mass and sea level](../docs/simulation-model/surface-processes.md#eustatic-sea-level-eustasypy)
- **Map projection lab.** Students compare the same coastline in Behrmann and Eckert IV,
  both equal-area projections, and explain how its apparent shape changes with the
  projection or viewing orientation.
  [Map projections](../docs/simulation-model/rendering-and-inspectors.md#projections-projectionspy)
- **Coastal planning hearing.** Teams recommend a site for a road or school under
  two modeled sea levels, then revisit the decision given information about access,
  cost, and residents. The teacher supplies the community constraints.
- **Computational modeling elective.** Students compare two model rules or parameter
  choices, then use a shared run record to explain the consequences. Synchronized maps
  and a small export of measurements could help students report reproducible comparisons.
