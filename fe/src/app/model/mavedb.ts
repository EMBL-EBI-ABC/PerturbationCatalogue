export interface MaveDBData {
  urn: string;
  title: string;
  shortDescription: string;
  sequenceType: string;
  geneName: string;
  geneCategory: string;
  publicationUrl: string;
  publicationYear: number;
  numVariants: number;
}

export interface MaveDBResponse {
  total: number;
  start: number;
  size: number;
  results: MaveDBData[];
  aggregations: any;
}

export interface MaveDBDetailsResponse {
  results: MaveDBData[];
}

export interface MaveDBFilters {
  sequenceType: string[];
  geneCategory: string[];
  publicationYear: number[];
}

