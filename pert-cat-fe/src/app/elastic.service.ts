import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { map } from 'rxjs/operators';


@Injectable({
  providedIn: 'root',
})
export class ElasticService {
  private readonly elasticUrl = 'https://USERNAME:KEY@SUBDOMAIN.elastic-cloud.com/mavedb/_search?filter_path=hits.hits._source';

  constructor(private http: HttpClient) { }

  getData(): Observable<any> {
    const query = {
      "size": 10,
      "query": {
        "match_all": {}
      }
    };

    const headers = new HttpHeaders().set('Content-Type', 'application/json');

    return this.http.post(this.elasticUrl, query, { headers }).pipe(
      map((response: any) => response.hits.hits.map((hit: { _source: any; }) => hit._source))
    );
  }
}
