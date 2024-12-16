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

  getData(start: number = 0, size: number = 10, query: string = ""): Observable<any> {
    const params = { start: start.toString(), size: size.toString(), q: query.toString() };
    return this.http.get(this.fastApiUrl, { params }).pipe(
      map((response: any) => ({
        results: response.results,
        total: response.total,
        aggregations: response.aggregations,
      }))
    );
  }
}
