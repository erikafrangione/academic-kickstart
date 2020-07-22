---
# Course title, summary, and position.
linktitle: Best Practices for Bioinformatics
summary: Sharing my thoughts on topics related to bioinformatics.
weight: 1

# Page metadata.
title: Data Management Best Practices for Informatics
date: "2018-09-09T00:00:00Z"
lastmod: "2018-09-09T00:00:00Z"
draft: false  # Is this a draft? true/false
toc: true  # Show table of contents? true/false
type: docs  # Do not modify.

# Add menu entry to sidebar.
# - name: Declare this menu item as a parent with ID `name`.
# - weight: Position of link in menu.
menu:
  example:
    name: Data Management
    weight: 1
---

If there’s one thing in my research that I struggled with as I began my journey into the realm of academia, it was data management. As a biology student, most of my prior experience with managing data involved hastily recording the amount of buffer I added into a chart from a ripped-out sheet of my lab manual, or plotting out a bacterial growth curve in Excel only for the original dataset to mysteriously vanish the day before my lab report was due. Most of my hands-on education up until grad school was heavily focused on learning procedures and interpreting results, and not so much on how to protect, store, and update the valuable data I was generating. I can’t count the number of times I used to throw my hands up in frustration over having lost or misplaced documents, or having a desktop flooded with filenames such as “dataset_v2.xlsx,” “dataset_final.xlsx,” “dataset_final_cleaned.xlsx,” or “dataset_final_edited.xlsx.” Never save important files to your desktop! On top of that, often these files would be hidden across multiple directories labeled “JUNK,” or “SCHOOL STUFF,” and sometimes leading to subdirectories with even more confusing and ambiguous filenames. I realized that there had to be a better way to manage my projects, especially when I was suddenly thrown into a role dealing with terabytes of sensitive and costly data. Good data management quickly became a priority within my research methods, and further, it allowed me to ensure my research is reliable and reproducible in the future.

Here, I’d like to share some of the best practices that I’ve learned for data management. Most of these come directly from my own experiences, and often failures, as well as from the helpful input of other researchers I’ve talked to. While there’s no directly standardized practice of managing data, I feel that I’ve come up with a reliable system that works for me as a researcher and aspiring bioinformatician. There are many different systems out there, but I feel these practices help me to stay organized and keep the focus on conducting good research.

## 1.	The Importance of Reproducible Research

Reproducible research should be a crucial component of any project or research undertaking. Any individual should be able to download my data, repeat my analyses, and still end up with the same results. This helps to ensure research credibility, strengthens the importance of our gathered evidence, extends the reach of our findings to other collaborators who may be interested in continuing or designing a related project, and allows for the ease of reusing existing scripts or project protocols.
One way we can ensure that our project is reproducible is by focusing on data management and organization from the get-go. Before even beginning our research, we can come up with a list of general questions to help us consider the types of data and data collection methods we might be dealing with, the different analyses we might need to perform, and the types of output files we might end up with. For example:
1.	How much raw data will we need to collect?
2.	Will all our raw data be accessible at once or collected over a period of time?
3.	How is our raw data structured and stored?
4.	How will we process our raw data (i.e. does it need to be cleaned, trimmed, manipulated, combined with other datasets, etc.)?
5.	What common analyses may we need to perform on our data?
6.	What are our desired outcomes for this project (i.e. are we looking to generate a specific data visualization, record the results from a specific analysis, add to our existing findings or interpret entirely new evidence, etc.)?

Developing a simple schematic for the project can help you in maintaining a structured and well-organized research management system. While many of your initial answers to these questions may change throughout the course of the research, understanding a project’s fundamental design can help you and fellow researchers ensure the results remain reproducible, and aid you in achieving your original project goals.

## Flexibility

This feature can be used for publishing content such as:

* **Online courses**
* **Project or software documentation**
* **Tutorials**

The `courses` folder may be renamed. For example, we can rename it to `docs` for software/project documentation or `tutorials` for creating an online course.

## Delete tutorials

**To remove these pages, delete the `courses` folder and see below to delete the associated menu link.**

## Update site menu

After renaming or deleting the `courses` folder, you may wish to update any `[[main]]` menu links to it by editing your menu configuration at `config/_default/menus.toml`.

For example, if you delete this folder, you can remove the following from your menu configuration:

```toml
[[main]]
  name = "Courses"
  url = "courses/"
  weight = 50
```

Or, if you are creating a software documentation site, you can rename the `courses` folder to `docs` and update the associated *Courses* menu configuration to:

```toml
[[main]]
  name = "Docs"
  url = "docs/"
  weight = 50
```

## Update the docs menu

If you use the *docs* layout, note that the name of the menu in the front matter should be in the form `[menu.X]` where `X` is the folder name. Hence, if you rename the `courses/example/` folder, you should also rename the menu definitions in the front matter of files within `courses/example/` from `[menu.example]` to `[menu.<NewFolderName>]`.

