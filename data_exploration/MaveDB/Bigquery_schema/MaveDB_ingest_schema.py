table_schema = {
    "fields": [
        {
            "name": "title",
            "type": "STRING",
            "mode": "REQUIRED"
        },
        {
            "name": "asOf",
            "type": "TIMESTAMP",
            "mode": "REQUIRED"
        },
        {
            "name": "experimentSets",
            "type": "RECORD",
            "mode": "REPEATED",
            "fields": [
                {
                    "name": "urn",
                    "type": "STRING",
                    "mode": "REQUIRED"
                },
                {
                    "name": "publishedDate",
                    "type": "DATE",
                    "mode": "REQUIRED"
                },
                {
                    "name": "id",
                    "type": "INTEGER",
                    "mode": "REQUIRED"
                },
                {
                    "name": "experiments",
                    "type": "RECORD",
                    "mode": "REPEATED",
                    "fields": [
                        {
                            "name": "title",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "shortDescription",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "abstractText",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "methodText",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "extraMetadata",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "keyname",
                                    "type": "STRING",
                                    "mode": ""
                                }
                            ]
                        },
                        {
                            "name": "urn",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "createdBy",
                            "type": "RECORD",
                            "mode": "REQUIRED",
                            "fields": [
                                {
                                    "name": "orcidId",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "firstName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "lastName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                }
                            ]
                        },
                        {
                            "name": "modifiedBy",
                            "type": "RECORD",
                            "mode": "REQUIRED",
                            "fields": [
                                {
                                    "name": "orcidId",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "firstName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "lastName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                }
                            ]
                        },
                        {
                            "name": "creationDate",
                            "type": "DATE",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "modificationDate",
                            "type": "DATE",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "publishedDate",
                            "type": "DATE",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "experimentSetUrn",
                            "type": "STRING",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "doiIdentifiers",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "identifier",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "id",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "url",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                }
                            ]
                        },
                        {
                            "name": "primaryPublicationIdentifiers",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "identifier",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "dbName",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "url",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "referenceHtml",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "title",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "abstract",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "authors",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "name",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "primary",
                                            "type": "BOOLEAN",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "publicationDoi",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "preprintDoi",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "publicationYear",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "preprintDate",
                                    "type": "DATE",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "publicationJournal",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "id",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                }
                            ]
                        },
                        {
                            "name": "secondaryPublicationIdentifiers",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "identifier",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "dbName",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "url",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "referenceHtml",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "title",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "abstract",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "authors",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "name",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "primary",
                                            "type": "BOOLEAN",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "publicationDoi",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "preprintDoi",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "publicationYear",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "preprintDate",
                                    "type": "DATE",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "publicationJournal",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "id",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                }
                            ]
                        },
                        {
                            "name": "rawReadIdentifiers",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "identifier",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "id",
                                    "type": "INTEGER",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "url",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                }
                            ]
                        },
                        {
                            "name": "keywords",
                            "type": "STRING",
                            "mode": "REPEATED"
                        },
                        {
                            "name": "scoreSets",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "title",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "methodText",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "abstractText",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "shortDescription",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "extraMetadata",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "dataUsagePolicy",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "urn",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "numVariants",
                                    "type": "INTEGER",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "license",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "longName",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "shortName",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "link",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "version",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "id",
                                            "type": "INTEGER",
                                            "mode": "REQUIRED"
                                        }
                                    ]
                                },
                                {
                                    "name": "supersededScoreSetUrn",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "supersedingScoreSetUrn",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                },
                                {
                                    "name": "metaAnalyzesScoreSetUrns",
                                    "type": "STRING",
                                    "mode": "REPEATED"
                                },
                                {
                                    "name": "metaAnalyzedByScoreSetUrns",
                                    "type": "STRING",
                                    "mode": "REPEATED"
                                },
                                {
                                    "name": "doiIdentifiers",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "identifier",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "id",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "url",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "primaryPublicationIdentifiers",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "identifier",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "dbName",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "url",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "referenceHtml",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "title",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "abstract",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "authors",
                                            "type": "RECORD",
                                            "mode": "REPEATED",
                                            "fields": [
                                                {
                                                    "name": "name",
                                                    "type": "STRING",
                                                    "mode": "NULLABLE"
                                                },
                                                {
                                                    "name": "primary",
                                                    "type": "BOOLEAN",
                                                    "mode": "NULLABLE"
                                                }
                                            ]
                                        },
                                        {
                                            "name": "publicationDoi",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "preprintDoi",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "publicationYear",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "preprintDate",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "publicationJournal",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "id",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "secondaryPublicationIdentifiers",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "identifier",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "dbName",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "url",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "referenceHtml",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "title",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "abstract",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "authors",
                                            "type": "RECORD",
                                            "mode": "REPEATED",
                                            "fields": [
                                                {
                                                    "name": "name",
                                                    "type": "STRING",
                                                    "mode": "NULLABLE"
                                                },
                                                {
                                                    "name": "primary",
                                                    "type": "BOOLEAN",
                                                    "mode": "NULLABLE"
                                                }
                                            ]
                                        },
                                        {
                                            "name": "publicationDoi",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "preprintDoi",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "publicationYear",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "preprintDate",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "publicationJournal",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "id",
                                            "type": "INTEGER",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "publishedDate",
                                    "type": "DATE",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "creationDate",
                                    "type": "DATE",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "modificationDate",
                                    "type": "DATE",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "createdBy",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "orcidId",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "firstName",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        },
                                        {
                                            "name": "lastName",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "modifiedBy",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "orcidId",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "firstName",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "lastName",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        }
                                    ]
                                },
                                {
                                    "name": "targetGenes",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "name",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "category",
                                            "type": "STRING",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "externalIdentifiers",
                                            "type": "RECORD",
                                            "mode": "REPEATED",
                                            "fields": [
                                                {
                                                    "name": "identifier",
                                                    "type": "RECORD",
                                                    "mode": "REPEATED",
                                                    "fields": [
                                                        {
                                                            "name": "dbName",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        },
                                                        {
                                                            "name": "identifier",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        },
                                                        {
                                                            "name": "dbVersion",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        },
                                                        {
                                                            "name": "url",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        },
                                                        {
                                                            "name": "referenceHtml",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        }
                                                    ]
                                                },
                                                {
                                                    "name": "offset",
                                                    "type": "INTEGER",
                                                    "mode": "REQUIRED"
                                                },
                                            ]
                                        },
                                        {
                                            "name": "id",
                                            "type": "INTEGER",
                                            "mode": "REQUIRED"
                                        },
                                        {
                                            "name": "targetSequence",
                                            "type": "RECORD",
                                            "mode": "REPEATED",
                                            "fields": [
                                                {
                                                    "name": "sequenceType",
                                                    "type": "STRING",
                                                    "mode": "REQUIRED"
                                                },
                                                {
                                                    "name": "sequence",
                                                    "type": "STRING",
                                                    "mode": "REQUIRED"
                                                },
                                                {
                                                    "name": "label",
                                                    "type": "STRING",
                                                    "mode": "NULLABLE"
                                                },
                                                {
                                                    "name": "taxonomy",
                                                    "type": "RECORD",
                                                    "mode": "REPEATED",
                                                    "fields": [
                                                        {
                                                            "name": "taxId",
                                                            "type": "INTEGER",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "organismName",
                                                            "type": "STRING",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "commonName",
                                                            "type": "STRING",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "rank",
                                                            "type": "STRING",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "hasDescribedSpeciesName",
                                                            "type": "BOOLEAN",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "articleReference",
                                                            "type": "STRING",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "genomeId",
                                                            "type": "STRING",
                                                            "mode": "NULLABLE"
                                                        },
                                                        {
                                                            "name": "id",
                                                            "type": "INTEGER",
                                                            "mode": "REQUIRED"
                                                        },
                                                        {
                                                            "name": "url",
                                                            "type": "STRING",
                                                            "mode": "REQUIRED"
                                                        }
                                                    ]
                                                }
                                            ]
                                        },
                                        {
                                            "name": "targetAccession",
                                            "type": "STRING",
                                            "mode": "NULLABLE"
                                        }
                                    ]
                                },
                                {
                                    "name": "datasetColumns",
                                    "type": "RECORD",
                                    "mode": "REPEATED",
                                    "fields": [
                                        {
                                            "name": "countColumns",
                                            "type": "STRING",
                                            "mode": "REPEATED"
                                        },
                                        {
                                            "name": "scoreColumns",
                                            "type": "STRING",
                                            "mode": "REPEATED"
                                        }
                                    ]
                                },
                                {
                                    "name": "keywords",
                                    "type": "STRING",
                                    "mode": "REPEATED"
                                },
                                {
                                    "name": "private",
                                    "type": "BOOLEAN",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "processingState",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "processingErrors",
                                    "type": "STRING",
                                    "mode": "NULLABLE"
                                }
                            ]
                        },
                        {
                            "name": "createdBy",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "orcidId",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "firstName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "lastName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                }
                            ]
                        },
                        {
                            "name": "modifiedBy",
                            "type": "RECORD",
                            "mode": "REPEATED",
                            "fields": [
                                {
                                    "name": "orcidId",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "firstName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                },
                                {
                                    "name": "lastName",
                                    "type": "STRING",
                                    "mode": "REQUIRED"
                                }
                            ]
                        },
                        {
                            "name": "creationDate",
                            "type": "DATE",
                            "mode": "REQUIRED"
                        },
                        {
                            "name": "modificationDate",
                            "type": "DATE",
                            "mode": "REQUIRED"
                        }
                    ]
                }
            ]
        }
    ]
}