import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';

@Injectable({
  providedIn: 'root',
})
export class ElasticService {
  private readonly fastApiUrl = 'https://perturbation-catalogue-be-959149465821.europe-west2.run.app/search';

  constructor(private http: HttpClient) {}

  getData(start: number = 0, size: number = 10, query: string = "", filters?: {
    publicationYear: any[];
    geneCategory: any[];
    sequenceType: any[]
  }): Observable<any> {
    let include_filters = false;
    // @ts-ignore
    let filter = '';
    // @ts-ignore
    if (filters.publicationYear.length > 0) {
      include_filters = true;
      // @ts-ignore
      filter = `publicationYear:${filters.publicationYear}`;
    }
    // @ts-ignore
    if (filters.geneCategory.length > 0) {
      if (include_filters) {
        // @ts-ignore
        filter = `${filter}+geneCategory:${filters.geneCategory}`;
      } else {
        include_filters = true;
        // @ts-ignore
        filter = `geneCategory:${filters.geneCategory}`;
      }
    }
    // @ts-ignore
    if (filters.sequenceType.length > 0) {
      if (include_filters) {
        // @ts-ignore
        filter = `${filter}+sequenceType:${filters.sequenceType}`;
      } else {
        include_filters = true;
        // @ts-ignore
        filter = `${filter}sequenceType:${filters.sequenceType}`;
      }
    }
    let params = { start: start.toString(), size: size.toString()};
    if (query) {
      // @ts-ignore
      params['q'] = query;
    }
    if (include_filters) {
      console.log(filter);
      // @ts-ignore
      params['filter'] = filter;
    }
    return this.http.get(this.fastApiUrl, { params }).pipe(
      map((response: any) => ({
        results: response.results,
        total: response.total,
        aggregations: response.aggregations,
      }))
    );
  }
}
