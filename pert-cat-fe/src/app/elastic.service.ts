import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';


@Injectable({
  providedIn: 'root',
})
export class ElasticService {
  private readonly fastApiUrl = 'https://pert-cat-be-959149465821.europe-west2.run.app/search';

  constructor(private http: HttpClient) { }

  getData(): Observable<any> {

    // Temporary substitution with first 10 entries from GitHub.
    const response = this.http.get(this.fastApiUrl)

    // Return the fields we need.
    return response.pipe(
      map((response: any) => response.results)
    );
  }
}
