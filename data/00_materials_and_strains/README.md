# Materials and strains

This folder contains the study-wide chemicals, materials, equipment, strains and medium formulations. Markdown tables and YAML files describe the materials and experimental preparations.

| Subject | Tables | Data |
|---|---|---|
| Chemicals, materials, and equipment used in this study | [chemicals_materials_equipment.md](chemicals_materials_equipment.md) | [chemicals_materials_equipment.yaml](chemicals_materials_equipment.yaml) |
| Strains used in this study | [strains.md](strains.md) | [strains.yaml](strains.yaml) |
| Media formulations used in this study | [media_formulations.md](media_formulations.md) | [media_formulations.yaml](media_formulations.yaml) |

The inventory contains 83 chemical/reagent-preparation entries (58 tested additives, 4 excluded preparations and 21 medium/assay reagents or preparations), 32 materials/equipment entries, 5 strains and 13 medium formulations.

This folder provides the shared inventories and definitions. Family YAMLs reference its recipes, strains, instruments and shared reference methods; experiment-specific assignments are documented within each family.

`SM4_BS_microplate` contains 0.19 g/L citric acid; `SM4_BS_tube` contains no citric acid. Culture assignments are defined in the experimental metadata.

The [stock preparation records](chemicals_materials_equipment.md#stock-preparation-records) contain reagent masses and final volumes for 14 preparations: 10 nominal additive stocks, the thymol blue stock and 3 excluded preparations. Preparation records are nested under the applicable `maximum_additive_stocks` or `attempted_stock_preparations` entries in [the inventory YAML](chemicals_materials_equipment.yaml), and under `reagent_preparations.thymol_blue` in [experimental_definitions.yaml](experimental_definitions.yaml).

Mass-derived concentrations use the formula mass without a purity correction. Preparation amounts and calculated nominal concentrations describe the nominal stocks. Excluded preparations list nominal preparation targets.

Contributor identities and roles are recorded in the [dataset overview](../README.md#contributors) and the relevant family READMEs. Return to the [data overview](../README.md).

## Reading shared YAML references

A mapping containing only `$ref` selects a value from another YAML file. The path is relative to the file containing the reference; the part after `#` is a JSON Pointer into that file. For example, `../00_materials_and_strains/strains.yaml#/strains/AF/strain` selects the AF strain designation. A pointer can select a scalar, list or full definition. Plain YAML parsers show the reference object; follow it to read the value.

[experimental_definitions.yaml](experimental_definitions.yaml) contains medium recipes, shared reference methods, reagent preparations and material identities. Local `medium_id`, `species_id`, physical-vessel, culture, acquisition and treatment assignments are defined in their experimental folders.

[chemical_reagents.yaml](chemical_reagents.yaml) exposes the `reagents` mapping through references to the bottle records in the materials inventory.

Keep this shared folder with the experimental data when transferring a subset of the deposit. Resolve each `$ref` by opening its relative YAML path and selecting the stated JSON Pointer.
