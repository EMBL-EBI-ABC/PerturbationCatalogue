import { Injectable } from '@angular/core';
import {HttpClient} from "@angular/common/http";
import {MaveDBDetailsResponse, MaveDBFilters, MaveDBResponse} from "../model/mavedb";
import {Observable} from "rxjs";

@Injectable({
  providedIn: 'root'
})
export class ApiService {

  constructor(private _httpClient: HttpClient) { }

  getAllMaveDBData(
    start: number,
    size: number,
    sortField: string = 'publicationYear',
    sortOrder: string = 'desc',
    query: any = null,
    filters: MaveDBFilters): Observable<MaveDBResponse> {
    const requestUrl = 'https://perturbation-catalogue-be-959149465821.europe-west2.run.app/search';
    let params = {
      start: start*size,
      size: size,
      sort_field: sortField,
      sort_order: sortOrder};

    if (query) {
      // @ts-ignore
      params['q'] = query;
    }

    let include_filters = false;
    let filter = '';
    if (filters.publicationYear.length > 0) {
      include_filters = true;
      filter = `publicationYear:${filters.publicationYear}`;
    }
    if (filters.geneCategory.length > 0) {
      if (include_filters) {
        filter = `${filter}+geneCategory:${filters.geneCategory}`;
      } else {
        include_filters = true;
        filter = `geneCategory:${filters.geneCategory}`;
      }
    }
    if (filters.sequenceType.length > 0) {
      if (include_filters) {
        filter = `${filter}+sequenceType:${filters.sequenceType}`;
      } else {
        include_filters = true;
        filter = `sequenceType:${filters.sequenceType}`;
      }
    }
    if (include_filters) {
      // @ts-ignore
      params['data_filter'] = filter;
    }
    return this._httpClient.get<MaveDBResponse>(requestUrl, {params: params});
  }

  getMaveDBRecord(recordId: string): Observable<MaveDBDetailsResponse> {
    const requestUrl = `https://perturbation-catalogue-be-959149465821.europe-west2.run.app/search/${recordId}`;
    return this._httpClient.get<MaveDBDetailsResponse>(requestUrl);
  }
}
